# -*- coding: utf-8 -*-
"""
MÓDULOS III E V: PURIFICACIÓN E DESTILACIÓN RIGOROSA (separation.py)
--------------------------------------------------------------------
Este módulo contén as unidades de separación:
1. DynamicFlash: Separador gas-líquido con nivel.
2. RecycleSplitter: Divisor de reciclo con lazo PID (mitigación Snowball).
3. DistillationColumn: Columnas simplificadas para Metanol (S2-1) e DME (D3-1, D3-2).
"""

import numpy as np
from control import PIDController

class DynamicFlash:
    """
    Separador flash con balance dinámico de líquido.
    """
    def __init__(self, T_set, P_set):
        self.T = T_set # K
        self.P = P_set # bar
        self.Level = 0.5 # 0 a 1
        
        # PID controlando o nivel a través do caudal de líquido
        self.pid_level = PIDController(Kp=1.0, Ki=0.1, Kd=0.0, Ts=2.0, action_type="direct", name="LIC-Flash")

    def compute_separation(self, feed_mol_h):
        """
        Calcula as correntes de vapor e líquido en equilibrio aproximado.
        feed_mol_h orde: [H2, CO, CO2, CH3OH, H2O, DME, N2, O2, CH4]
        """
        # Simplificación de flash: gases permanentes van ao vapor, líquidos ao líquido
        # segundo un k-value empírico baseado en volatilidades á T e P dadas.
        # Vapor pressure aproximado (bar) a T=323 K (50 C)
        # H2, CO, N2, CH4, O2 son moi volátiles.
        k_values = np.array([1000.0, 1000.0, 100.0, 0.15, 0.05, 5.0, 1000.0, 1000.0, 1000.0])
        
        # Clamp feed to be non-negative and check for NaN or zero-sum
        feed_mol_h = np.maximum(feed_mol_h, 0.0)
        total_feed = np.sum(feed_mol_h)
        if total_feed < 1e-8 or np.isnan(total_feed):
            return np.zeros_like(feed_mol_h), np.zeros_like(feed_mol_h)
            
        z = feed_mol_h / total_feed
        V_F = 0.5
        for _ in range(10): # Newton-Raphson rápido
            denom = 1.0 + V_F * (k_values - 1.0)
            denom = np.where(np.abs(denom) < 1e-10, 1e-10 * np.sign(denom), denom)
            f = np.sum(z * (1.0 - k_values) / denom)
            df = np.sum(-z * (1.0 - k_values)**2 / denom**2)
            if np.abs(df) < 1e-12:
                break
            V_F = V_F - f / df
            V_F = np.clip(V_F, 0.001, 0.999)
            
        V_total = V_F * total_feed
        L_total = (1.0 - V_F) * total_feed
        
        x = z / (1.0 + V_F * (k_values - 1.0))
        y = k_values * x
        
        vapor_mol_h = y * V_total
        liquid_mol_h = x * L_total
        
        return vapor_mol_h, liquid_mol_h

    def step_level(self, liquid_in_mol_h, dt=2.0):
        """
        Actualiza o nivel e devolve o factor de apertura de saída de líquido.
        """
        u_valve_liq = self.pid_level.compute(sp=0.5, pv=self.Level)
        
        # Balance volumétrico aproximado
        flow_in_m3_s = np.sum(liquid_in_mol_h) * 0.03 / 3600.0 # aprox 30 cm3/mol
        flow_out_m3_s = u_valve_liq * 0.1 # max 0.1 m3/s
        
        dV_dt = flow_in_m3_s - flow_out_m3_s
        self.Level += dV_dt * dt / 10.0 # 10 m3 de tanque
        self.Level = np.clip(self.Level, 0.0, 1.0)
        
        return u_valve_liq


class RecycleSplitter:
    """
    SPT2-1: Divisor de reciclo con control PID.
    """
    def __init__(self, target_split=0.75):
        self.target_split = target_split
        # O PID controla o purga para evitar acumulación excesiva (Snowball)
        # Se a presión do sistema ou inventario sobe moito, aumenta a purga.
        self.pid_purge = PIDController(Kp=0.5, Ki=0.05, Kd=0.0, Ts=2.0, action_type="direct", name="PIC-Splitter")
        self.split_fraction = target_split
        # Novo PID de composición para Nitróxeno
        self.pid_n2 = PIDController(Kp=1.0, Ki=0.05, Kd=0.0, Ts=2.0, action_type="direct", u_min=-0.15, u_max=0.15, name="AIC-N2")

    def split(self, feed_mol_h, system_pressure=10.0):
        """
        Recircula unha fracción ao reactor, purga o resto. (Lóxica baseada en presión)
        """
        purge_valve = self.pid_purge.compute(sp=10.0, pv=system_pressure)
        purge_fraction = (1.0 - self.target_split) + (purge_valve - 0.5) * 0.2
        purge_fraction = np.clip(purge_fraction, 0.01, 0.99)
        self.split_fraction = 1.0 - purge_fraction
        
        recycle = feed_mol_h * self.split_fraction
        purge = feed_mol_h * purge_fraction
        
        return recycle, purge

    def split_dynamic_purge(self, feed_mol_h, n2_in_mol_h):
        """
        Calcula a purga baseada no balance estrito de N2 para estabilizar en 18.54% (Mitigación do Snowball Effect)
        """
        total_feed = np.sum(feed_mol_h)
        if total_feed > 0:
            y_N2 = feed_mol_h[6] / total_feed
        else:
            y_N2 = 0.0
            
        # Fracción base necesaria para purgar os moles exactos que entraron con este paso temporal
        if feed_mol_h[6] > 1e-5:
            f_purge_base = n2_in_mol_h / feed_mol_h[6]
        else:
            f_purge_base = 0.25
            
        # Controlador para asegurar a calibración arredor do set point (0.1854)
        purge_adj = self.pid_n2.compute(sp=0.1854, pv=y_N2)
        
        purge_fraction = f_purge_base + purge_adj
        purge_fraction = np.clip(purge_fraction, 0.001, 0.999)
        self.split_fraction = 1.0 - purge_fraction
        
        recycle = feed_mol_h * self.split_fraction
        purge = feed_mol_h * purge_fraction
        
        return recycle, purge, y_N2


class DistillationColumn:
    """
    Columna de destilación con modelo simplificado (Masa e Enerxía global).
    Permite establecer calidade por cabeza e controlar o calderín.
    """
    def __init__(self, n_stages, feed_stage, RR, Q_reboiler, Q_condenser):
        self.N = n_stages
        self.F_stage = feed_stage
        self.RR = RR
        self.Q_reboiler = Q_reboiler
        self.Q_condenser = Q_condenser
        
        # PID de calidade por cabeza (Control de pureza axustando Refluxo ou Destilado)
        self.pid_top = PIDController(Kp=5.0, Ki=0.1, Kd=0.0, Ts=2.0, action_type="reverse", name="TC-Top")
        # PID de fondos (Control do calderín)
        self.pid_bot = PIDController(Kp=10.0, Ki=0.5, Kd=0.0, Ts=2.0, action_type="direct", name="TC-Bot")
        
        self.T_top = 40.0
        self.T_bot = 150.0

    def compute(self, feed_mol_h, purity_target_top, T_set_bot):
        """
        Realiza un balance MESH simplificado baseado na recuperación de chaves (Key Components).
        """
        # Simplificando a complexidade do solver ríxido, usamos un modelo de separación 
        # con axuste de compoñentes lixeiros á cabeza e pesados a colas, influenciado pola enerxía.
        
        feed_mol_h = np.maximum(feed_mol_h, 0.0)
        total_feed = np.sum(feed_mol_h)
        if total_feed < 1e-5 or np.isnan(total_feed):
            return np.zeros_like(feed_mol_h), np.zeros_like(feed_mol_h)
            
        z = feed_mol_h / total_feed
        
        # Compoñentes: [H2, CO, CO2, CH3OH, H2O, DME, N2, O2, CH4]
        # Volatilidades relativas aproximadas: H2>N2>CO>O2>CH4>CO2>DME>CH3OH>H2O
        alpha = np.array([100.0, 80.0, 50.0, 0.8, 0.1, 2.5, 90.0, 70.0, 60.0])
        
        # Controlador axusta un parámetro de "sharpness" (nitidez da separación)
        sharpness_top = self.pid_top.compute(sp=purity_target_top, pv=self.T_top) # Usamos T_top como proxy temporal ou asumimos medidor directo
        Q_reboiler_adj = self.pid_bot.compute(sp=T_set_bot, pv=self.T_bot)
        
        # Recuperación á cabeza de cada compoñente
        # En vez de función contínua, aplicamos un corte case perfecto para compoñentes lixeiros vs pesados.
        recovery_top = np.where(alpha >= 1.0, 0.999, 0.001)
        
        top_mol_h = feed_mol_h * recovery_top
        bot_mol_h = feed_mol_h * (1.0 - recovery_top)
        
        # Actualizar temperaturas estimadas
        y_top = top_mol_h / np.maximum(np.sum(top_mol_h), 1e-5)
        x_bot = bot_mol_h / np.maximum(np.sum(bot_mol_h), 1e-5)
        
        # Estimación de T de burbulla (moi groseira)
        # Baseada nas propiedades dos maioritarios
        self.T_top = np.sum(y_top * np.array([-250, -190, -78, 64, 100, -24, -195, -183, -161])) + 273.15 # T_eb C -> K
        self.T_bot = np.sum(x_bot * np.array([-250, -190, -78, 64, 100, -24, -195, -183, -161])) + 273.15
        
        return top_mol_h, bot_mol_h
