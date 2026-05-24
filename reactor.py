# -*- coding: utf-8 -*-
"""
MÓDULO CINÉTICO E DE REACTOR DINÁMICO (reactor.py)
--------------------------------------------------
Este módulo modela o reactor catalítico de leito fixo para a síntese directa de DME.
- Cinética LHHW para síntese de Metanol, WGS e deshidratación de Metanol a DME.
- Balances de materia e enerxía discretizados (1D por método de liñas).
- Modelado de transferencia térmica co fluído de refrixeración (coolant).
"""

import numpy as np
from thermodynamics import R, TC, PC, OMEGA, PengRobinsonEoS

# Propiedades térmicas aproximadas (J/mol*K)
CP_GAS = np.array([29.0, 29.0, 37.0, 44.0, 33.0, 65.0]) # [H2, CO, CO2, CH3OH, H2O, DME]
CP_CAT = 800.0   # J/(kg_cat * K) - calor específico do catalizador sólido

# Entalpías de reacción a 298K (J/mol)
# R1: CO + 2H2 <-> CH3OH         (Exotérmica)
# R2: CO + H2O <-> CO2 + H2      (Exotérmica, WGS)
# R3: 2 CH3OH <-> DME + H2O      (Exotérmica, Deshidratación)
DH_RXN = np.array([-90.8e3, -41.2e3, -23.5e3])

# Coeficientes estequiométricos (nu_ij)
# Mistura: [H2, CO, CO2, CH3OH, H2O, DME]
# Reacción 1 (MeOH): -2*H2 - 1*CO + 1*CH3OH
# Reacción 2 (WGS):  +1*H2 - 1*CO + 1*CO2 - 1*H2O
# Reacción 3 (DME):  -2*CH3OH + 1*H2O + 1*DME
NU = np.array([
    [-2.0, -1.0,  0.0,  1.0,  0.0,  0.0],  # R1
    [ 1.0, -1.0,  1.0,  0.0, -1.0,  0.0],  # R2
    [ 0.0,  0.0,  0.0, -2.0,  1.0,  1.0]   # R3
])

class LHHWKinetics:
    """
    Cinética de Langmuir-Hinshelwood-Hougen-Watson (LHHW) para a síntese directa de DME.
    As presións parciais deben estar en bar.
    """
    def __init__(self):
        # Constantes pre-exponenciais e enerxías de activación (Ea en J/mol)
        # k = k0 * exp(-Ea / (R*T))
        # R1: Síntese de Metanol
        self.k1_0 = 1.2e4      # mol/(kg_cat * s * bar^2)
        self.Ea1 = 65.0e3      # J/mol
        
        # R2: Water-Gas Shift (WGS)
        self.k2_0 = 8.5e3      # mol/(kg_cat * s * bar)
        self.Ea2 = 55.0e3      # J/mol
        
        # R3: Deshidratación a DME
        self.k3_0 = 3.5e5      # mol/(kg_cat * s * bar^2)
        self.Ea3 = 80.0e3      # J/mol

        # Constantes de adsorción (K_i = K0_i * exp(-DH_ads_i / (R*T)))
        self.K_CO_0 = 0.5      # bar^-1
        self.DH_CO = -15e3     # J/mol
        
        self.K_H2_0 = 0.1      # bar^-0.5
        self.DH_H2 = -10e3     # J/mol
        
        self.K_MeOH_0 = 1.2    # bar^-1
        self.DH_MeOH = -25e3   # J/mol

        self.K_H2O_0 = 2.0     # bar^-1
        self.DH_H2O = -30e3    # J/mol

    def compute_rates(self, T, P_bar, y, phi):
        """
        Calcula as velocidades de reacción (r1, r2, r3) en mol/(kg_cat * s).
        Utiliza as presións parciais corrixidas por fugacidade (fuxacidades): f_i = y_i * P * phi_i
        """
        # Presións parciais en bar (fuxacidades para ser termodinamicamente consistentes)
        f = y * P_bar * phi
        
        f = np.clip(f, 1e-10, 1000.0) # Evitar valores negativos
        f_H2, f_CO, f_CO2, f_MeOH, f_H2O, f_DME = f

        # Constantes de equilibrio químico K_eq = exp(A + B/T) (P en bar)
        # R1: CO + 2H2 <-> CH3OH
        K_eq1 = np.exp(-21.22 + 10540.0 / T)
        # R2: CO + H2O <-> CO2 + H2
        K_eq2 = np.exp(-4.33 + 4577.8 / T)
        # R3: 2 CH3OH <-> DME + H2O
        K_eq3 = np.exp(-9.56 + 3383.0 / T)

        # Calcular constantes cinéticas
        k1 = self.k1_0 * np.exp(-self.Ea1 / (R * T))
        k2 = self.k2_0 * np.exp(-self.Ea2 / (R * T))
        k3 = self.k3_0 * np.exp(-self.Ea3 / (R * T))

        # Calcular constantes de adsorción
        K_CO = self.K_CO_0 * np.exp(-self.DH_CO / (R * T))
        K_H2 = self.K_H2_0 * np.exp(-self.DH_H2 / (R * T))
        K_MeOH = self.K_MeOH_0 * np.exp(-self.DH_MeOH / (R * T))
        K_H2O = self.K_H2O_0 * np.exp(-self.DH_H2O / (R * T))

        # Reacción 1: Síntese de Metanol
        # Num = f_CO * f_H2^2 - f_MeOH / K_eq1
        num1 = f_CO * (f_H2**2) - (f_MeOH / K_eq1)
        den1 = (1.0 + K_CO * f_CO + (K_H2 * f_H2)**0.5 + K_MeOH * f_MeOH)**3
        r1 = k1 * num1 / den1

        # Reacción 2: WGS
        # Num = f_CO * f_H2O - f_CO2 * f_H2 / K_eq2
        num2 = f_CO * f_H2O - (f_CO2 * f_H2 / K_eq2)
        den2 = (1.0 + K_CO * f_CO + K_H2O * f_H2O + 0.5 * f_CO2)**2
        r2 = k2 * num2 / den2

        # Reacción 3: Deshidratación a DME
        # Num = f_MeOH^2 - f_DME * f_H2O / K_eq3
        num3 = f_MeOH**2 - (f_DME * f_H2O / K_eq3)
        den3 = (1.0 + K_MeOH * f_MeOH + K_H2O * f_H2O)**2
        r3 = k3 * num3 / den3

        return np.array([r1, r2, r3])


class FixedBedReactor:
    """
    Representa un reactor de leito fixo catalítico discretizado en N_nodes.
    O estado do reactor está definido por:
      - n_i,k: moles de cada compoñente i no nodo k (N_nodes * 6 variables)
      - T_k: temperatura no nodo k (N_nodes variables)
    Total de variables de estado: N_nodes * 7
    """
    def __init__(self, N_nodes=10, L=3.0, D_inner=0.1, mass_cat_total=500.0):
        self.N = N_nodes
        self.L = L
        self.D = D_inner
        self.Area = np.pi * (D_inner**2) / 4.0
        self.V_total = self.Area * L
        self.V_cell = self.V_total / N_nodes
        self.dz = L / N_nodes
        
        self.mass_cat_total = mass_cat_total
        self.W_cat_cell = mass_cat_total / N_nodes
        
        self.kinetics = LHHWKinetics()
        self.eos = PengRobinsonEoS()
        
        # Parámetros de transporte
        self.U_heat = 150.0  # Coeficiente global de transferencia de calor (W / m2*K)
        self.A_heat_cell = np.pi * D_inner * self.dz  # Área de transferencia por cell (m2)
        
        # Parámetros de dinámica hidráulica
        # Resistencia ao fluxo entre celdas (molar flow = K_drag * deltaP)
        self.K_drag = 1.5e-4  # mol/(s * Pa)
        
        # Presión de saída por defecto (bar)
        self.P_outlet_set = 50.0

    def compute_derivatives(self, t, state, F_feed, T_feed, T_coolant, u_valve, P_sink):
        """
        Calcula as derivadas dState/dt.
        state: array unidimensional de tamaño N*7
          [ n_0,0, n_1,0, ... n_5,0, T_0,    n_0,1, ... n_5,1, T_1,   ... ]
        F_feed: caudal molar de alimentación de compoñentes [mol/s] (lonxitude 6)
        T_feed: temperatura de alimentación [K]
        T_coolant: temperatura do refrixerante [K] (pode ser un escalar ou un array de tamaño N)
        u_valve: apertura da válvula de saída (0 a 1)
        P_sink: presión augas abaixo do reactor (Pa)
        """
        N = self.N
        W_cell = self.W_cat_cell
        
        # Desempaquetar o estado
        # n[k, i] = moles do compoñente i no nodo k
        # T[k] = temperatura no nodo k
        n = np.zeros((N, 6))
        T = np.zeros(N)
        for k in range(N):
            n[k, :] = state[k*7 : k*7 + 6]
            T[k] = state[k*7 + 6]
            
        # Evitar moles negativos por oscilacións numéricas
        n = np.clip(n, 1e-6, 1e6)
        T = np.clip(T, 200.0, 1000.0)
        
        # Calcular variables auxiliares para cada nodo: total moles, presións, fugacidades, etc.
        n_total = np.sum(n, axis=1)
        y = np.zeros((N, 6))
        P = np.zeros(N)
        phi = np.zeros((N, 6))
        Z = np.zeros(N)
        
        for k in range(N):
            y[k, :] = n[k, :] / n_total[k]
            # Usar PR-EoS para estimar compresibilidade Z e coeficientes de fugacidade phi
            # Primeiro facemos unha estimación ideal da presión
            P_est = (n_total[k] * R * T[k]) / self.V_cell
            Z_k, phi_k = self.eos.compute_mixture_properties(T[k], P_est, y[k, :])
            Z[k] = Z_k
            phi[k, :] = phi_k
            # Presión real corrixida polo factor Z
            P[k] = (n_total[k] * Z_k * R * T[k]) / self.V_cell
            
        # Determinar os fluxos molares interactivos F_k entre os nodos
        # F[k, i] = fluxo molar do compoñente i que sae do nodo k e entra no nodo k+1
        F = np.zeros((N, 6))
        
        # Fluxo de entrada ao primeiro nodo (k=0) é o feed
        F_in = F_feed # mol/s
        
        # Fluxos entre nodos k -> k+1 (para k = 0 a N-2)
        for k in range(N - 1):
            dP = P[k] - P[k+1]
            F_total_k = self.K_drag * dP  # mol/s
            
            if F_total_k >= 0:
                # Fluxo directo (cara a adiante)
                F[k, :] = y[k, :] * F_total_k
            else:
                # Fluxo inverso (cara a atrás)
                F[k, :] = y[k+1, :] * F_total_k
                
        # Fluxo de saída do último nodo (k = N-1) controlado pola válvula de control de presión
        # F_out = C_valve * u_valve * sqrt(P_last - P_sink)
        C_valve = 8e-3  # Constante de capacidade da válvula (mol / (s * Pa^0.5))
        dP_valve = np.max([0.0, P[N-1] - P_sink])
        F_total_valve = C_valve * u_valve * np.sqrt(dP_valve)
        F[N-1, :] = y[N-1, :] * F_total_valve
        
        # Calcular velocidades de reacción e derivadas de estado para cada nodo
        dn_dt = np.zeros((N, 6))
        dT_dt = np.zeros(N)
        
        # Temperatura do refrixerante por nodo
        if np.isscalar(T_coolant):
            Tc_array = np.full(N, T_coolant)
        else:
            Tc_array = T_coolant
            
        for k in range(N):
            # Fluxo de entrada ao nodo k
            F_node_in = F_feed if k == 0 else F[k-1, :]
            # Fluxo de saída do nodo k
            F_node_out = F[k, :]
            
            # Presión en bar
            P_bar = P[k] / 1e5
            
            # Velocidades de reacción
            rates = self.kinetics.compute_rates(T[k], P_bar, y[k, :], phi[k, :])
            
            # Consumo/produción química de compoñentes
            # r_i = sum_j (nu_ij * r_j)
            rxn_term = np.zeros(6)
            for i in range(6):
                rxn_term[i] = np.sum(NU[:, i] * rates)
                
            # Derivada de materia: dn_i/dt = F_in_i - F_out_i + W_cat * rxn_term_i
            dn_dt[k, :] = F_node_in - F_node_out + W_cell * rxn_term
            
            # Enerxía: Inercia térmica integrada (catalizador + gas)
            # cp_mix_gas = sum(y_i * cp_gas_i)
            cp_gas_mix = np.sum(y[k, :] * CP_GAS)
            thermal_mass = n_total[k] * cp_gas_mix + W_cell * CP_CAT
            
            # Calor xerado polas reaccións (exotérmicas)
            Q_rxn = W_cell * np.sum(-DH_RXN * rates)
            
            # Calor transferido co refrixerante
            Q_cool = self.U_heat * self.A_heat_cell * (T[k] - Tc_array[k])
            
            # Ence de enerxía convectivo
            # sum(F_in_i * cp_i * (T_in - T_k))
            T_in = T_feed if k == 0 else T[k-1]
            Q_conv = np.sum(F_node_in * CP_GAS * (T_in - T[k]))
            
            # dT/dt = (Q_conv + Q_rxn - Q_cool) / thermal_mass
            dT_dt[k] = (Q_conv + Q_rxn - Q_cool) / thermal_mass
            
        # Aplanar derivadas para retornar
        derivatives = np.zeros(N * 7)
        for k in range(N):
            derivatives[k*7 : k*7 + 6] = dn_dt[k, :]
            derivatives[k*7 + 6] = dT_dt[k]
            
        return derivatives


# Proba rápida de funcionamento do reactor
if __name__ == "__main__":
    reactor = FixedBedReactor(N_nodes=5)
    
    # Alimentación molar típica (mol/s): H2: 6.5, CO: 2.5, CO2: 0.5, CH3OH: 0.0, H2O: 0.0, DME: 0.0
    F_feed = np.array([6.5, 2.5, 0.5, 0.0, 0.0, 0.0]) # Total 9.5 mol/s
    T_feed = 523.15  # 250 C
    T_coolant = 513.15  # 240 C
    u_valve = 0.5
    P_sink = 10.0 * 1e5  # 10 bar
    
    # Inicializar estado uniforme (Z=1 aproximado para establecer moles iniciais a 50 bar)
    P_init = 50.0 * 1e5
    n_init_node = (F_feed / np.sum(F_feed)) * (P_init * reactor.V_cell) / (R * T_feed)
    state = np.zeros(5 * 7)
    for k in range(5):
        state[k*7 : k*7 + 6] = n_init_node
        state[k*7 + 6] = T_feed
        
    derivs = reactor.compute_derivatives(0.0, state, F_feed, T_feed, T_coolant, u_valve, P_sink)
    print("TEST REACTOR:")
    print(f"Número de variables de estado: {len(state)}")
    print(f"Derivadas calculadas no primeiro nodo (dn_dt): {derivs[0:6]}")
    print(f"Derivada de temperatura no primeiro nodo (dT_dt): {derivs[6]:.4f} K/s")
