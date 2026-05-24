# -*- coding: utf-8 -*-
"""
MÓDULO II-A: COMPRESIÓN SEGURA MULTIETAPA (compression.py)
----------------------------------------------------------
Este módulo simula o compresor con límite mecánico de relación por etapa e intercooling.
Inclúe o quecedor previo H2-1.
"""

import numpy as np
from control import PIDController

class MultiStageCompressor:
    """
    Compresor Isentrópico Multietapa (K2-1) con quentador previo (H2-1).
    Restrición mecánica: A relación de compresión por etapa non pode superar 3.
    """
    def __init__(self):
        # Quentador previo
        self.T_in_H2_1 = 40.0 + 273.15 # K
        self.T_out_H2_1 = 120.0 + 273.15 # K
        
        # Compresor
        self.P_in = 1.0 # bar
        self.P_out_setpoint = 110.0 # bar
        self.T_intercool = 150.0 + 273.15 # K
        self.ratio_max = 3.0
        
        # Rendemento isentrópico
        self.eta_is = 0.85
        
        # Controlador PID de velocidade do motor (RPM) para manter a presión de descarga
        self.pid_speed = PIDController(Kp=1.5, Ki=0.2, Kd=0.1, Ts=2.0, 
                                       u_min=500.0, u_max=3000.0, action_type="reverse", name="SIC-K2")
        
        # O controlador reverse: se P_out cae por baixo do setpoint, sp-pv > 0 -> aumenta u (velocidade).
        # Agarda: se P_out cae, queremos aumentar velocidade para presurizar máis, entón o erro é SP - PV.
        # Direct: PV - SP. Se PV > SP (sobrepresión), baixa RPM.
        # Polo tanto action_type debe ser "reverse" (SP - PV)
        
        self.rpm = 1500.0 # RPM inicial
        self.P_out_current = self.P_in

    def calculate_stages(self, P_in, P_out):
        """
        Calcula o número de etapas e a relación de compresión para non superar ratio_max=3
        """
        if P_out <= P_in:
            return 1, 1.0
            
        ratio_total = P_out / P_in
        # ratio_total = ratio_stage ^ n_stages => n_stages = ln(ratio_total) / ln(ratio_stage)
        # Para que ratio_stage <= 3: n_stages >= ln(ratio_total) / ln(3)
        n_stages_min = np.ceil(np.log(ratio_total) / np.log(self.ratio_max))
        n_stages = int(max(1, n_stages_min))
        
        ratio_stage = ratio_total ** (1.0 / n_stages)
        return n_stages, ratio_stage

    def compute(self, syngas_mol_h, P_in, T_in, target_P_out=110.0):
        """
        Calcula as temperaturas de descarga reais e o traballo consumido usando leis isentrópicas
        T_out_isentropic = T_in * (P_out/P_in) ^ ((gamma-1)/gamma)
        """
        # H2-1: Quentar o gas
        T_gas = self.T_out_H2_1
        
        # Propiedades térmicas medias simplificadas do syngas (H2, CO)
        gamma = 1.4 # Coeficiente adiabático aproximado para diatómicos
        
        n_stages, ratio_stage = self.calculate_stages(P_in, target_P_out)
        
        work_total_W = 0.0
        T_current = T_gas
        P_current = P_in
        
        flow_mol_s = np.sum(syngas_mol_h) / 3600.0
        R = 8.314 # J/molK
        Cp_avg = R * gamma / (gamma - 1.0)
        
        for i in range(n_stages):
            P_next = P_current * ratio_stage
            
            # Temperatura de descarga isentrópica
            T_out_is = T_current * (ratio_stage ** ((gamma - 1.0)/gamma))
            
            # Temperatura real baseada na eficiencia
            T_out_real = T_current + (T_out_is - T_current) / self.eta_is
            
            # Traballo da etapa
            W_stage = flow_mol_s * Cp_avg * (T_out_real - T_current)
            work_total_W += W_stage
            
            # Intercooling entre etapas (agás na última, aínda que ás veces si)
            if i < n_stages - 1:
                T_current = self.T_intercool
            else:
                T_current = T_out_real # T_descarga
                
            P_current = P_next
            
        # O controlador de velocidade (RPM) modula lixeiramente o fluxo ou a relación real 
        # en dinámica, pero neste modelo usaremos o PID para fixar o RPM que mantería a presión.
        # Simulamos que a presión de saída é controlada por este RPM.
        
        return P_current, T_current, work_total_W, n_stages, ratio_stage

    def step_dynamics(self, pv_P_out, dt=2.0):
        """
        Calcula o novo valor de RPM en función da presión lida.
        """
        self.P_out_current = pv_P_out
        self.rpm = self.pid_speed.compute(sp=self.P_out_setpoint, pv=self.P_out_current)
        return self.rpm
