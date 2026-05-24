# -*- coding: utf-8 -*-
"""
MÓDULO II-B: REACTOR DE SÍNTESE DE METANOL (methanol_reactor.py)
----------------------------------------------------------------
Este módulo modela o reactor PFR multitubular de leito fixo (R2-1)
con cinética LHHW heteroxénea estrita para as reaccións R18 (Síntese) e R19 (WGS Inversa).
"""

import numpy as np
from thermodynamics import R, PengRobinsonEoS
from control import PIDController

class MethanolLHHWKinetics:
    """
    Cinética LHHW para Síntese de Metanol (R18) e WGS Inversa (R19).
    Velocidades en kmol/(kg_cat * s).
    """
    def __init__(self):
        pass

    def compute_rates(self, T, P_bar, y, phi):
        """
        T: K, P_bar: bar, y: fraccións molares, phi: coeficientes de fugacidade.
        y orde: [H2, CO, CO2, CH3OH, H2O, DME, N2, O2, CH4]
        """
        # Presións parciais reais (fugacidades) en bar
        f = y * P_bar * phi
        f = np.clip(f, 1e-10, 1000.0)
        
        P_H2 = f[0]
        P_CO = f[1]
        P_CO2 = f[2]
        P_MeOH = f[3]
        P_H2O = f[4]
        
        # Constantes de velocidade k1 e k2 para R18 e R19
        k1_18 = 1.07e-13 * np.exp(4413.76 / T)
        k2_18 = 4.182e7 * np.exp(-2645.97 / T)
        
        k1_19 = 122.0 * np.exp(-11398.24 / T)
        k2_19 = 1.1412 * np.exp(-6624.98 / T)
        
        # Constantes de adsorción K1, K2, K3
        K1 = 3453.38
        K2 = 1.578e-3 * np.exp(2068.44 / T)
        K3 = 6.62e-16 * np.exp(14928.92 / T)
        
        # Denominador común
        den_base = 1.0 + K1 * (P_H2O / P_H2) + K2 * np.sqrt(P_H2) + K3 * P_H2O
        
        # R18: CO2 + 3H2 <-> CH3OH + H2O
        num_18 = k1_18 * P_CO * P_H2 - k2_18 * (P_MeOH * P_H2O / np.sqrt(P_H2))
        r18 = num_18 / den_base
        
        # R19: CO2 + H2 <-> CO + H2O
        num_19 = k1_19 * P_CO2 - k2_19 * (P_CO * P_H2O / np.sqrt(P_H2))
        r19 = num_19 / np.sqrt(den_base)
        
        return np.array([r18, r19])


class MethanolPFR:
    """
    Reactor PFR de leito fixo 1D.
    Resolve masa e enerxía discretizando espacialmente (Método das Liñas).
    """
    def __init__(self, N_nodes=10):
        self.N = N_nodes
        
        # Xeometría R2-1
        self.N_tubes = 2000
        self.L = 0.35 # m
        self.D_tube = 0.015 # m
        self.eps = 0.5 # void fraction
        self.rho_cat = 2000.0 # kg/m3 (partícula)
        self.rho_bed = 1000.0 # kg/m3 (aparente)
        
        # Volumes e Áreas
        self.Area_tube = np.pi * (self.D_tube**2) / 4.0
        self.V_tube = self.Area_tube * self.L
        self.V_total = self.V_tube * self.N_tubes
        
        self.V_cell = self.V_total / self.N
        self.W_cat_total = 123.7 # kg
        self.W_cell = self.W_cat_total / self.N
        self.dz = self.L / self.N
        
        self.kinetics = MethanolLHHWKinetics()
        self.eos = PengRobinsonEoS()
        
        # Entalpías (J/kmol)
        self.DH_R18 = -41000.0 * 1000.0
        self.DH_R19 = 51000.0 * 1000.0
        
        # Capacidades caloríficas gas (J/kmol*K) aprox
        self.CP_GAS = np.array([29.0, 29.0, 37.0, 44.0, 33.0, 65.0, 29.0, 29.0, 35.0]) * 1000.0
        self.CP_CAT = 800.0 # J/kg*K
        
        # Transferencia térmica
        self.U_heat = 150.0 # W/m2*K
        self.A_heat_cell = np.pi * self.D_tube * self.dz * self.N_tubes
        
        # Estequiometría: [H2, CO, CO2, CH3OH, H2O, DME, N2, O2, CH4]
        # R18: -3 H2, -1 CO2, +1 CH3OH, +1 H2O
        # R19: -1 H2, -1 CO2, +1 CO, +1 H2O
        self.NU = np.array([
            [-3.0,  0.0, -1.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0], # R18
            [-1.0,  1.0, -1.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0]  # R19
        ])
        
        # Controladores PID de nivel de planta
        self.pid_temp = PIDController(Kp=1.0, Ki=0.01, Kd=0.1, Ts=2.0, action_type="reverse", name="TC-R2")
        self.pid_temp = PIDController(Kp=1.0, Ki=0.01, Kd=0.1, Ts=2.0, action_type="reverse", name="TC-R2")
        self.T_coolant = 200.0 + 273.15 # K
        
        self.K_drag = 2.0e-7 # kmol/(s * Pa) - Coeficiente para fluxo

    def compute_derivatives(self, t, state, F_feed_kmols, T_feed, u_valve, P_sink):
        """
        state: [n_i_0, T_0, n_i_1, T_1 ... ] n en kmol, T en K. 10 variables por nodo (9 comp + 1 T)
        """
        N = self.N
        n = np.zeros((N, 9))
        T = np.zeros(N)
        
        for k in range(N):
            n[k, :] = state[k*10 : k*10 + 9]
            T[k] = state[k*10 + 9]
            
        n_pos = np.maximum(n, 0.0)
        T = np.clip(T, 200.0, 1000.0)
        
        n_total = np.sum(n_pos, axis=1)
        y = np.zeros((N, 9))
        P = np.zeros(N)
        phi = np.zeros((N, 9))
        Z = np.zeros(N)
        
        for k in range(N):
            n_tot_k = max(n_total[k], 1e-10)
            y[k, :] = n_pos[k, :] / n_tot_k
            P_est = (n_tot_k * 1000.0 * R * T[k]) / self.V_cell # P = nRT/V en Pa (n en mol)
            Z_k, phi_k = self.eos.compute_mixture_properties(T[k], P_est, y[k, :])
            Z[k] = Z_k
            phi[k, :] = phi_k
            P[k] = (n_tot_k * 1000.0 * Z_k * R * T[k]) / self.V_cell
            
        F = np.zeros((N, 9))
        for k in range(N - 1):
            dP = P[k] - P[k+1]
            F_total_k = self.K_drag * dP # kmol/s
            if F_total_k >= 0:
                F[k, :] = y[k, :] * F_total_k
            else:
                F[k, :] = y[k+1, :] * F_total_k
                
        # Válvula final
        C_valve = 5e-7
        dP_valve = np.max([0.0, P[N-1] - P_sink])
        F_total_valve = C_valve * u_valve * np.sqrt(dP_valve)
        F[N-1, :] = y[N-1, :] * F_total_valve
        
        dn_dt = np.zeros((N, 9))
        dT_dt = np.zeros(N)
        
        for k in range(N):
            F_in = F_feed_kmols if k == 0 else F[k-1, :]
            F_out = F[k, :]
            
            P_bar = P[k] / 1e5
            rates = self.kinetics.compute_rates(T[k], P_bar, y[k, :], phi[k, :])
            
            rxn_term = np.zeros(9)
            for i in range(9):
                rxn_term[i] = np.sum(self.NU[:, i] * rates)
                
            dn_dt[k, :] = F_in - F_out + self.W_cell * rxn_term
            
            # Enerxía
            cp_gas_mix = np.sum(y[k, :] * self.CP_GAS)
            thermal_mass = n_total[k] * cp_gas_mix + self.W_cell * self.CP_CAT
            
            Q_rxn = self.W_cell * (-self.DH_R18 * rates[0] - self.DH_R19 * rates[1])
            Q_cool = self.U_heat * self.A_heat_cell * (T[k] - self.T_coolant)
            
            T_in = T_feed if k == 0 else T[k-1]
            Q_conv = np.sum(F_in * self.CP_GAS * (T_in - T[k]))
            
            dT_dt[k] = (Q_conv + Q_rxn - Q_cool) / thermal_mass
            
        derivatives = np.zeros(N * 10)
        for k in range(N):
            derivatives[k*10 : k*10 + 9] = dn_dt[k, :]
            derivatives[k*10 + 9] = dT_dt[k]
            
        return derivatives
