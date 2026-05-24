# -*- coding: utf-8 -*-
"""
MÓDULO TERMODINÁMICO E DE PROPIEDADES FÍSICAS (thermodynamics.py)
------------------------------------------------------------------
Este módulo contén a lóxica para os cálculos termodinámicos e propiedades físicas.
Fase Gas: Ecuación de Estado de Peng-Robinson (PR-EoS) para misturas de H2, CO, CO2, CH3OH, H2O, DME.
Fase Líquida: Modelo de coeficientes de actividade NRTL para DME, Metanol e Auga.
"""

import numpy as np

# R = 8.314462618 J/(mol*K) - Constante universal dos gases ideais
R = 8.314462618

# Propiedades críticas dos compoñentes da mistura:
# Orde: [H2, CO, CO2, CH3OH, H2O, DME, N2, O2, CH4]
NAMES = ["H2", "CO", "CO2", "CH3OH", "H2O", "DME", "N2", "O2", "CH4"]
TC = np.array([33.19, 132.85, 304.13, 512.6, 647.1, 400.0, 126.2, 154.6, 190.6])       # K
PC = np.array([13.13, 34.94, 73.77, 80.97, 220.64, 52.4, 34.00, 50.43, 45.99]) * 1e5    # Pa (bar * 1e5)
OMEGA = np.array([-0.216, 0.048, 0.224, 0.564, 0.344, 0.200, 0.037, 0.022, 0.011])      # Factor acéntrico

# Matriz de parámetros de interacción binaria (kij) para Peng-Robinson (estimacións típicas)
# H2, CO, CO2, CH3OH, H2O, DME, N2, O2, CH4
KIJ_PR = np.array([
    [0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],   # H2
    [0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],   # CO
    [0.0,   0.0,   0.0,   0.05,  0.05,  0.05,  0.0,   0.0,   0.09],  # CO2
    [0.0,   0.0,   0.05,  0.0,   -0.07, 0.02,  0.0,   0.0,   0.0],   # CH3OH
    [0.0,   0.0,   0.05,  -0.07, 0.0,   0.08,  0.0,   0.0,   0.0],   # H2O
    [0.0,   0.0,   0.05,  0.02,  0.08,  0.0,   0.0,   0.0,   0.0],   # DME
    [0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.03],  # N2
    [0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0],   # O2
    [0.0,   0.0,   0.09,  0.0,   0.0,   0.0,   0.03,  0.0,   0.0]    # CH4
])

class PengRobinsonEoS:
    """
    Ecuación de Estado de Peng-Robinson para misturas multi-compoñente.
    """
    def __init__(self, Tc=TC, Pc=PC, omega=OMEGA, kij=KIJ_PR):
        self.Tc = Tc
        self.Pc = Pc
        self.omega = omega
        self.kij = kij
        self.num_components = len(Tc)
        
        # Parámetro kappa para cada compoñente
        # kappa_i = 0.37464 + 1.54226*omega_i - 0.26992*omega_i^2
        self.kappa = 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2
        
        # Parámetros a_c e b para Peng-Robinson
        self.a_c = 0.45724 * (R**2 * self.Tc**2) / self.Pc
        self.b = 0.07780 * (R * self.Tc) / self.Pc

    def compute_mixture_properties(self, T, P, y):
        """
        Calcula o factor de compresibilidade Z e os coeficientes de fugacidade (phi)
        para unha mestura de gas a temperatura T (K), presión P (Pa) e fraccións molares y.
        """
        # Evitar valores de temperatura e presión non físicos ou nulos
        T = max(T, 100.0)
        P = max(P, 1e3)
        
        # Evitar división por cero ou fraccións negativas en y
        if np.any(np.isnan(y)) or np.sum(y) < 1e-8:
            y = np.zeros(self.num_components)
            y[6] = 1.0 # 100% N2 por defecto
        else:
            y = np.maximum(y, 0.0)
            y = y / np.sum(y)
        
        # 1. Calcular a_i(T) para cada compoñente
        Tr = T / self.Tc
        alpha = (1.0 + self.kappa * (1.0 - np.sqrt(Tr)))**2
        a = self.a_c * alpha
        
        # 2. Regras de mestura (van der Waals)
        # b_mix = sum(y_i * b_i)
        b_mix = np.sum(y * self.b)
        b_mix = max(b_mix, 1e-6)
        
        # a_mix = sum_i sum_j (y_i * y_j * sqrt(a_i * a_j) * (1 - kij))
        a_matrix = np.zeros((self.num_components, self.num_components))
        for i in range(self.num_components):
            for j in range(self.num_components):
                a_matrix[i, j] = np.sqrt(a[i] * a[j]) * (1.0 - self.kij[i, j])
                
        a_mix = 0.0
        for i in range(self.num_components):
            a_mix += y[i] * np.sum(y * a_matrix[i, :])
        a_mix = max(a_mix, 1e-6)
            
        # 3. Parámetros adimensionais da ecuación de estado
        A = (a_mix * P) / (R**2 * T**2)
        B = (b_mix * P) / (R * T)
        
        # 4. Resolver a ecuación cúbica en Z: Z^3 - (1 - B)*Z^2 + (A - 2B - 3B^2)*Z - (AB - B^2 - B^3) = 0
        coeffs = [
            1.0,
            -(1.0 - B),
            A - 2.0*B - 3.0*B**2,
            -(A*B - B**2 - B**3)
        ]
        
        if np.any(np.isnan(coeffs)) or np.any(np.isinf(coeffs)):
            return 1.0, np.ones(self.num_components)
            
        try:
            roots = np.roots(coeffs)
            real_roots = roots[np.isreal(roots)].real
            # Filtrar valores físicos (Z > B)
            physical_roots = real_roots[real_roots > B]
            if len(physical_roots) > 0:
                Z = np.max(physical_roots)
            else:
                Z = 1.0
        except Exception:
            Z = 1.0
            
        # 5. Cálculo de coeficientes de fugacidade (phi)
        # ln(phi_i) = b_i/b * (Z - 1) - ln(Z - B) - A / (2*sqrt(2)*B) * (2*sum_j(y_j * a_ij)/a - b_i/b) * ln((Z + (1+sqrt(2))B)/(Z + (1-sqrt(2))B))
        phi = np.zeros(self.num_components)
        
        # Termo logarítmico común
        denom_log = Z + (1.0 - np.sqrt(2.0)) * B
        if denom_log <= 0:
            denom_log = 1e-10
        log_term = np.log(max(1e-10, Z + (1.0 + np.sqrt(2.0)) * B) / denom_log)
        
        log_Z_minus_B = np.log(max(1e-10, Z - B))
        
        for i in range(self.num_components):
            term_sum = np.sum(y * a_matrix[i, :])
            term_ratio = 2.0 * term_sum / a_mix
            
            ln_phi_i = (self.b[i] / b_mix) * (Z - 1.0) - log_Z_minus_B - (A / (2.0 * np.sqrt(2.0) * B)) * (term_ratio - self.b[i] / b_mix) * log_term
            # Evitar desbordamento de exp co limitador
            phi[i] = np.exp(np.clip(ln_phi_i, -11.5, 11.5))
            
        # Limitar phis para evitar problemas numéricos
        phi = np.clip(phi, 1e-5, 1e5)
        
        return Z, phi


class NRTLModel:
    """
    Modelo NRTL (Non-Random Two-Liquid) para misturas altamente non ideais.
    Usado para a fase líquida ternaria: DME (1) - Metanol (2) - Auga (3).
    """
    def __init__(self):
        # Parámetros NRTL: A_ij (K) = (g_ij - g_jj) / R
        # Fila i, Columna j correspondente a A_ij
        # Índices: 0 = DME, 1 = Metanol, 2 = Auga
        # Valores de literatura aproximados para destilación de DME:
        self.A = np.array([
            [0.0,      485.6,   1250.2],   # DME (0) -> [DME, MeOH, H2O]
            [-120.4,   0.0,     385.1],    # MeOH (1) -> [DME, MeOH, H2O]
            [850.1,    420.3,   0.0]       # H2O (2) -> [DME, MeOH, H2O]
        ])
        
        # Parámetros de non-aleatoriedade alpha_ij (alfa_ij = alfa_ji)
        self.alpha = np.array([
            [0.0,   0.3,   0.3],
            [0.3,   0.0,   0.3],
            [0.3,   0.3,   0.0]
        ])

    def compute_activity_coefficients(self, T, x):
        """
        Calcula os coeficientes de actividade (gamma) a unha temperatura T (K) e fraccións molares líquidas x.
        x debe ser de lonxitude 3: [DME, Metanol, Auga]
        """
        x = np.clip(x, 1e-15, 1.0)
        x = x / np.sum(x)
        
        num_comp = len(x)
        tau = self.A / T
        G = np.exp(-self.alpha * tau)
        
        gamma = np.zeros(num_comp)
        for i in range(num_comp):
            # Termo 1: sum_j(tau_ji * G_ji * x_j) / sum_k(G_ki * x_k)
            num1 = np.sum(tau[:, i] * G[:, i] * x)
            den1 = np.sum(G[:, i] * x)
            term1 = num1 / den1
            
            # Termo 2: sum_j [ (x_j * G_ij / sum_k(G_kj * x_k)) * (tau_ij - (sum_m(x_m * tau_mj * G_mj) / sum_k(G_kj * x_k))) ]
            term2 = 0.0
            for j in range(num_comp):
                den_j = np.sum(G[:, j] * x)
                term_sum_m = np.sum(x * tau[:, j] * G[:, j])
                
                term_inner = tau[i, j] - (term_sum_m / den_j)
                term2 += (x[j] * G[i, j] / den_j) * term_inner
                
            gamma[i] = np.exp(term1 + term2)
            
        return gamma


# Proba rápida de funcionamento
if __name__ == "__main__":
    # Test Peng-Robinson
    pr = PengRobinsonEoS()
    # Mistura molar típica de entrada ao reactor: H2: 0.65, CO: 0.25, CO2: 0.05, CH3OH: 0.01, H2O: 0.0, DME: 0.04, N2: 0.0, O2: 0.0, CH4: 0.0
    y_test = np.array([0.65, 0.25, 0.05, 0.01, 0.0, 0.04, 0.0, 0.0, 0.0])
    y_test = y_test / np.sum(y_test)
    T_test = 523.15  # 250 C
    P_test = 50.0 * 1e5  # 50 bar
    
    Z, phi = pr.compute_mixture_properties(T_test, P_test, y_test)
    print("TEST PENG-ROBINSON:")
    print(f"Compressibility Factor Z: {Z:.4f}")
    for name, ph in zip(NAMES, phi):
        print(f"  Fugacity coeff phi_{name}: {ph:.4f}")
        
    # Test NRTL
    nrtl = NRTLModel()
    x_test = np.array([0.2, 0.5, 0.3])  # DME, MeOH, H2O
    T_liq = 348.15  # 75 C
    gamma = nrtl.compute_activity_coefficients(T_liq, x_test)
    print("\nTEST NRTL:")
    for name, gam in zip(["DME", "MeOH", "H2O"], gamma):
        print(f"  Activity coeff gamma_{name}: {gam:.4f}")
