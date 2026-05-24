# -*- coding: utf-8 -*-
"""
MÓDULO I: FRONT-END, GASIFICACIÓN E LAZOS ATMOSFÉRICOS (gasification.py)
-----------------------------------------------------------------------
Este módulo contén os modelos para a sección de gasificación da biomasa:
1. DryerRStoic: Secadoiro dinámico.
2. PyrolysisRYield: Reactor de pirólise.
3. GasifierRGibbs: Reactor de gasificación baseado na minimización da Enerxía Libre de Gibbs.
4. Tren de arrefriamento e FlashSeparator con PID.
"""

import numpy as np
from scipy.optimize import minimize
from thermodynamics import R
from control import PIDController

class FeedstockManager:
    """
    Libraría dinámica de materiais de biomasa para a alimentación.
    """
    def __init__(self):
        self.materials = {
            "pine_woodchips": {"C": 0.4345, "H": 0.0424, "O": 0.3910, "N": 0.0072, "S": 0.0040, "Ash": 0.1209},
            "vine_shoots":    {"C": 0.4710, "H": 0.0580, "O": 0.4200, "N": 0.0090, "S": 0.0020, "Ash": 0.0400},
            "straw":          {"C": 0.4100, "H": 0.0500, "O": 0.3800, "N": 0.0100, "S": 0.0050, "Ash": 0.1450},
            "sawdust":        {"C": 0.4900, "H": 0.0600, "O": 0.4300, "N": 0.0050, "S": 0.0010, "Ash": 0.0140}
        }
        self.current_material = "pine_woodchips"
        self.custom_comp = None
        
    def set_material(self, name):
        if name in self.materials:
            self.current_material = name
            self.custom_comp = None
            
    def set_custom_composition(self, C, H, O, N, S, Ash):
        self.custom_comp = {"C": C, "H": H, "O": O, "N": N, "S": S, "Ash": Ash}
        self.current_material = "custom"
        
    def get_composition(self):
        if self.custom_comp:
            return self.custom_comp
        return self.materials[self.current_material]

# Pesos atómicos (kg/kmol)
MW_C = 12.011
MW_H = 1.008
MW_O = 15.999
MW_N = 14.007
MW_S = 32.065
MW_H2O = 18.015

class DryerRStoic:
    """
    Secadoiro Dinámico R1-1 (RStoic).
    Operación: 150 C, 1 bar.
    Converte a biomasa húmida liberando auga.
    """
    def __init__(self):
        self.T = 150.0 + 273.15 # K
        self.P = 1.0 # bar
        self.conversion = 0.117

    def compute(self, mass_flow_biomass):
        """
        mass_flow_biomass: kg/h de biomasa entrante
        Retorna: (kg/h biomasa seca acondicionada, kg/h auga liberada)
        """
        # Reacción: Biomasa -> 0.0555084 H2O
        # A estequiometría indica que o 11.7% da masa é humidade liberada
        water_released = mass_flow_biomass * self.conversion
        biomass_conditioned = mass_flow_biomass - water_released
        return biomass_conditioned, water_released


class PyrolysisRYield:
    """
    Reactor de Pirólise R1-2 (RYield).
    Operación: 750 C, 1 bar.
    Descompón a biomasa nos seus elementos e separa as cinzas.
    """
    def __init__(self):
        self.T = 750.0 + 273.15 # K
        self.P = 1.0 # bar

    def compute(self, mass_flow_biomass_dry, water_vapor_in, comp):
        """
        Calcula os fluxos molares atómicos resultantes.
        Retorna:
        - moles/h de [C, H, O, N, S]
        - kg/h de cinzas (sólido)
        - kg/h de vapor de auga que pasa intacto
        """
        ash_flow = mass_flow_biomass_dry * comp["Ash"]
        
        # O resto é masa orgánica que se descompón
        organic_mass = mass_flow_biomass_dry - ash_flow
        
        # Normalizar fraccións orgánicas
        sum_org = comp["C"] + comp["H"] + comp["O"] + comp["N"] + comp["S"]
        fC = comp["C"] / sum_org
        fH = comp["H"] / sum_org
        fO = comp["O"] / sum_org
        fN = comp["N"] / sum_org
        fS = comp["S"] / sum_org
        
        # Moles de átomos (kmol/h) * 1000 = mol/h
        mol_C = (organic_mass * fC / MW_C) * 1000.0
        mol_H = (organic_mass * fH / MW_H) * 1000.0
        mol_O = (organic_mass * fO / MW_O) * 1000.0
        mol_N = (organic_mass * fN / MW_N) * 1000.0
        mol_S = (organic_mass * fS / MW_S) * 1000.0
        
        # Engadir os átomos do vapor de auga entrante á corrente elemental
        mol_H2O_in = (water_vapor_in / MW_H2O) * 1000.0
        mol_H += 2.0 * mol_H2O_in
        mol_O += mol_H2O_in
        
        atoms_mol_h = np.array([mol_C, mol_H, mol_O, mol_N, mol_S])
        
        return atoms_mol_h, ash_flow, water_vapor_in


class GasifierRGibbs:
    """
    Gasificador de Equilibrio R1-3 (RGibbs) con Lazo de Temperatura.
    Minimiza Enerxía Libre de Gibbs.
    """
    def __init__(self):
        self.T_setpoint = 750.0 + 273.15 # K
        self.P = 1.0 # bar
        # PID para manter a temperatura controlando o fluxo de calor Q
        self.pid_temp = PIDController(Kp=50.0, Ki=2.0, Kd=5.0, Ts=2.0, 
                                      u_min=-1e7, u_max=1e7, action_type="reverse", name="TC-RGibbs")
        self.Q_heater = 0.0
        self.T_current = self.T_setpoint
        self.C_solid_out_kgh = 0.0

    def compute_equilibrium(self, atoms_mol_h, air_flow_kgh):
        """
        atoms_mol_h: [C, H, O, N, S] mol/h
        air_flow_kgh: kg/h (79% N2, 21% O2 en masa aprox. ou molar, asumimos molar para simplificar e axustamos masa)
        """
        # Converter aire (asumindo fracción molar 79% N2, 21% O2) a moles/h
        MW_air = 0.79 * (MW_N*2) + 0.21 * (MW_O*2)
        air_mol_h = (air_flow_kgh / MW_air) * 1000.0
        
        mol_N2_air = 0.79 * air_mol_h
        mol_O2_air = 0.21 * air_mol_h
        
        # Engadir átomos do aire
        total_C = atoms_mol_h[0]
        total_H = atoms_mol_h[1]
        total_O = atoms_mol_h[2] + mol_O2_air * 2.0
        total_N = atoms_mol_h[3] + mol_N2_air * 2.0
        
        # Gases de saída posibles: H2, CO, CO2, H2O, CH4, N2
        # Moles_saida: [n_H2, n_CO, n_CO2, n_H2O, n_CH4, n_N2]
        
        # Para evitar unha optimización lenta e non lineal complexa en cada paso, 
        # resolvemos un sistema simplificado de balances de masa e constantes de equilibrio 
        # (WGS e Metanación) que se achega ao RGibbs.
        # N2 é inerte.
        n_N2 = total_N / 2.0
        
        # Para garantir o balance atómico estrito, usamos só 2 graos de liberdade
        # Variables independentes: n_CH4, n_CO
        # Resto variables:
        # n_CO2 = total_C - n_CO - n_CH4
        # n_H2O = total_O - n_CO - 2.0*n_CO2
        # n_H2 = (total_H - 2.0*n_H2O - 4.0*n_CH4) / 2.0
        
        def equations(vars):
            n_CH4, n_CO = vars
            
            n_CO2 = total_C - n_CO - n_CH4
            n_H2O = total_O - n_CO - 2.0 * n_CO2
            n_H2 = (total_H - 2.0 * n_H2O - 4.0 * n_CH4) / 2.0
            
            all_vars = [n_H2, n_CO, n_CO2, n_H2O, n_CH4]
            
            if any(v < 0 for v in all_vars):
                return 1e9 + sum(abs(v) for v in all_vars if v < 0) * 1e6
                
            n_total = sum(all_vars) + n_N2
            
            T = self.T_current
            K_WGS = np.exp(4220.0/T - 3.86)
            K_Meth = np.exp(26830.0/T - 30.11)
            
            P_H2 = n_H2 / n_total
            P_CO = n_CO / n_total
            P_CO2 = n_CO2 / n_total
            P_H2O = n_H2O / n_total
            P_CH4 = n_CH4 / n_total
            
            err_eq1 = (P_CO2 * P_H2) - K_WGS * (P_CO * P_H2O)
            err_eq2 = (P_CH4 * P_H2O) - K_Meth * (P_CO * (P_H2**3))
            
            return err_eq1**2 * n_total**2 + err_eq2**2 * n_total**4
            
        x0 = [total_C*0.05, min(total_C*0.6, total_O*0.8)]
        bounds = [(0, total_H/4.0), (0, total_C)]
        
        res = minimize(equations, x0, method='SLSQP', bounds=bounds, options={'maxiter': 1000, 'ftol': 1e-8})
        
        n_CH4, n_CO = res.x
        n_CO2 = total_C - n_CO - n_CH4
        n_H2O = total_O - n_CO - 2.0 * n_CO2
        n_H2 = (total_H - 2.0 * n_H2O - 4.0 * n_CH4) / 2.0
        
        n_CH4 = max(0, n_CH4)
        n_CO = max(0, n_CO)
        n_CO2 = max(0, n_CO2)
        n_H2O = max(0, n_H2O)
        n_H2 = max(0, n_H2)
        
        # O sólido só o calculamos analíticamente despois para a gráfica
        self.C_solid_out_kgh = 0.0
        
        # Fluxo total en mol/h
        syngas_molar_flow = np.array([n_H2, n_CO, n_CO2, 0.0, n_H2O, 0.0, n_N2, 0.0, n_CH4])
        
        return syngas_molar_flow

    def step_thermal(self, syngas_flow_mol_h, dt=2.0):
        """
        Simula o lazo térmico
        """
        # Controlador de Temperatura
        self.Q_heater = self.pid_temp.compute(sp=self.T_setpoint, pv=self.T_current)
        # Dinámica de 1ª orde térmica simulada
        thermal_mass = 50000.0 # J/K aprox
        dT = (self.Q_heater) / thermal_mass * dt
        self.T_current += dT
        return self.T_current


class FlashSeparator:
    """
    Separación, Arrefriamento e Control Flash (SEP1-1, E1-1, E1-2, F1-1).
    F1-1 opera a 40 C e 1 bar.
    """
    def __init__(self):
        self.T_flash = 40.0 + 273.15 # K
        self.P_flash = 1.0 # bar
        # PID para presión usando válvula de gas
        self.pid_pressure = PIDController(Kp=0.1, Ki=0.01, Kd=0.05, Ts=2.0, action_type="direct", name="PIC-F1")
        # PID para nivel usando válvula de líquido
        self.pid_level = PIDController(Kp=2.0, Ki=0.1, Kd=0.0, Ts=2.0, action_type="direct", name="LIC-F1")
        
        self.P_current = 1.0
        self.Level_current = 0.5 # 0 a 1
        
    def compute(self, syngas_in_mol_h, ash_in_kgh):
        """
        syngas_in_mol_h: [H2, CO, CO2, CH3OH, H2O, DME, N2, O2, CH4]
        Realiza a separación de auga a 40C (condensación).
        A 40C, a presión de vapor da auga é baixa, polo que moita condensa.
        """
        # Impurezas sólidas (SEP1-1)
        ash_out = ash_in_kgh 
        
        # Arrefriamento (E1-1, E1-2 asumidos ideais que baixan T a 40C)
        
        # F1-1 Flash a 40C. 
        # Vapor pressure de H2O a 40C aprox 0.0738 bar (Antoine)
        Pv_H2O = 0.0738
        
        total_gas_dry = syngas_in_mol_h[0] + syngas_in_mol_h[1] + syngas_in_mol_h[2] + syngas_in_mol_h[6] + syngas_in_mol_h[8]
        
        # Fracción molar de auga máxima no gas
        y_H2O_max = Pv_H2O / self.P_flash
        
        # Moles de auga que saturan o gas:
        # n_H2O_gas / (total_gas_dry + n_H2O_gas) = y_H2O_max
        # n_H2O_gas = y_H2O_max * total_gas_dry / (1 - y_H2O_max)
        n_H2O_gas = y_H2O_max * total_gas_dry / (1.0 - y_H2O_max)
        
        n_H2O_total = syngas_in_mol_h[4]
        
        if n_H2O_total > n_H2O_gas:
            n_H2O_liquid = n_H2O_total - n_H2O_gas
            n_H2O_gas_out = n_H2O_gas
        else:
            n_H2O_liquid = 0.0
            n_H2O_gas_out = n_H2O_total
            
        syngas_out_mol_h = np.copy(syngas_in_mol_h)
        syngas_out_mol_h[4] = n_H2O_gas_out
        
        liquid_out_mol_h = np.zeros(9)
        liquid_out_mol_h[4] = n_H2O_liquid
        
        return syngas_out_mol_h, liquid_out_mol_h, ash_out
        
    def step_dynamics(self, syngas_out_mol_h, liquid_out_mol_h, dt=2.0):
        # O PIC regula a válvula de gas (saída 4) para manter 1 bar
        u_valve_gas = self.pid_pressure.compute(sp=1.0, pv=self.P_current)
        
        # O LIC regula a válvula de líquido (RLIQ) para manter nivel
        u_valve_liq = self.pid_level.compute(sp=0.5, pv=self.Level_current)
        
        # Ecuación de volume diferencial (simplificada)
        # Nivel sube con líquido entrante, baixa con válvula
        liq_flow_m3_h = liquid_out_mol_h[4] * 18.015 / 1000.0 / 1000.0 # m3/h aprox
        dV_dt = liq_flow_m3_h - (u_valve_liq * 0.5) 
        self.Level_current += dV_dt * (dt/3600.0) / 2.0 # Asumindo área 2m2
        
        return u_valve_gas, u_valve_liq
