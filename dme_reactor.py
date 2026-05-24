# -*- coding: utf-8 -*-
"""
MÓDULO IV: SECCIÓN DE SÍNTESE DE DME (dme_reactor.py)
-----------------------------------------------------
Este módulo contén a bomba P3-1, os lazos térmicos reversibles (H3-1, H3-2)
e o reactor de equilibrio de DME (R3-1).
"""

import numpy as np
from control import PIDController

class HydraulicPump:
    """
    Bomba P3-1: Eleva presión de líquido.
    """
    def __init__(self):
        self.P_in = 2.6 # bar
        self.P_out = 15.1 # bar
        
    def compute(self, flow_mol_h, P_in, P_out_set):
        # Asumimos incompresible, o traballo é mínimo (W = V * dP)
        # Non incluímos cálculo de traballo detallado por simplicidade.
        return P_out_set, flow_mol_h


class ReversibleHeater:
    """
    Lazo Térmico Reversible (H3-1, H3-2).
    A acción do PID está programada para que, se T > SP, peche a válvula 
    (acción inversa respecto a un arrefriador normal).
    Como o fluído quente (utilidade) quenta a corrente de proceso:
    T_process < SP -> Necesitamos máis calor -> Abrir válvula (aumentar u).
    Se T_process sube, o erro (SP - PV) baixa ou faise negativo, reducindo u.
    Isto é un PID "reverse acting" estándar onde u regula o caudal quente.
    Se fose un enfriador: T_process > SP -> Necesitamos máis frío -> Abrir válvula (direct acting).
    """
    def __init__(self, target_T, name="HC-3"):
        self.T_setpoint = target_T + 273.15
        self.pid_temp = PIDController(Kp=2.0, Ki=0.1, Kd=0.05, Ts=2.0, action_type="reverse", name=name)
        self.u_valve = 0.5
        self.T_current = self.T_setpoint - 10.0 # Empeza un pouco máis frío

    def step_thermal(self, T_in, dt=2.0):
        """
        Calcula a nova temperatura baseada na apertura da válvula de servizo.
        """
        self.u_valve = self.pid_temp.compute(sp=self.T_setpoint, pv=self.T_current)
        
        # O calor achegado pola utilidade é proporcional á válvula
        Q_utility = self.u_valve * 50000.0 # Watts máximos
        
        # Asumimos unha inercia térmica do equipo
        thermal_mass = 20000.0 # J/K
        dT_dt = (Q_utility - 500.0 * (self.T_current - T_in)) / thermal_mass
        self.T_current += dT_dt * dt
        
        return self.T_current, self.u_valve


class DMEEquilibriumReactor:
    """
    Reactor de Equilibrio de DME (R3-1).
    2 CH3OH <-> DME + H2O
    """
    def __init__(self):
        self.T = 250.0 + 273.15 # K
        self.P = 14.7 # bar
        
    def compute_equilibrium(self, feed_mol_h):
        """
        Calcula o equilibrio a 250 C calibrado segundo o target do obxectivo.
        feed_mol_h: [H2, CO, CO2, CH3OH, H2O, DME, N2, O2, CH4]
        """
        # Clamping and copying the feed to be robust against negative/NaN flows
        feed_mol_h = np.maximum(np.nan_to_num(feed_mol_h), 0.0)
        out_mol_h = np.copy(feed_mol_h)
        
        n_MeOH_in = feed_mol_h[3]
        n_H2O_in = feed_mol_h[4]
        n_DME_in = feed_mol_h[5]
        
        # Asumimos constante de equilibrio a 250C para a reacción: 2 MeOH <-> DME + H2O
        # Keq_250 = (P_DME * P_H2O) / P_MeOH^2
        # Aproximado para dar conversións típicas de ~80-85%
        Keq = 5.0 
        
        # Grao de avance (xi)
        # MeOH = n_MeOH_in - 2*xi
        # DME = n_DME_in + xi
        # H2O = n_H2O_in + xi
        
        # (n_DME_in + xi)*(n_H2O_in + xi) / (n_MeOH_in - 2*xi)^2 = Keq
        
        # Resolver a ecuación cuadrática para xi
        # (D+x)(H+x) = K*(M-2x)^2
        # DH + Dx + Hx + x^2 = K*(M^2 - 4Mx + 4x^2)
        # x^2 + (D+H)x + DH = 4Kx^2 - 4KMx + KM^2
        # (1 - 4K)x^2 + (D + H + 4KM)x + (DH - KM^2) = 0
        
        D = n_DME_in
        H = n_H2O_in
        M = n_MeOH_in
        
        a = 1.0 - 4.0 * Keq
        b = D + H + 4.0 * Keq * M
        c = D * H - Keq * M**2
        
        # Se a = 0 (Keq=0.25)
        if abs(a) < 1e-5:
            xi = -c / b
        else:
            discriminant = b**2 - 4*a*c
            if discriminant >= 0:
                xi1 = (-b + np.sqrt(discriminant)) / (2*a)
                xi2 = (-b - np.sqrt(discriminant)) / (2*a)
                # O xi válido é aquel que non dá moles negativos
                xi = xi1 if (M - 2*xi1 >= 0 and xi1 >= 0) else xi2
            else:
                xi = 0.0 # Non hai reacción
                
        # Protexer contra conversión excesiva
        max_xi = M / 2.0
        xi = np.clip(xi, 0.0, max_xi * 0.99)
        
        out_mol_h[3] = M - 2.0 * xi
        out_mol_h[4] = H + xi
        out_mol_h[5] = D + xi
        
        return out_mol_h
