# -*- coding: utf-8 -*-
"""
MÓDULO VI: MOTOR CENTRAL E INTEGRACIÓN PLANTWIDE (engine.py)
------------------------------------------------------------
Este módulo implementa a clase `DigitalTwinEngine`, que conecta 
todos os equipos e resolve o lazo de simulación (Simulation Loop).
Acopla os lazos de control e as integracións numéricas no tempo.
"""

import numpy as np
from scipy.integrate import solve_ivp

# Importar bloques modulares
from gasification import DryerRStoic, PyrolysisRYield, GasifierRGibbs, FlashSeparator, FeedstockManager
from gas_compression import MultiStageCompressor
from methanol_reactor import MethanolPFR
from separation import DynamicFlash, RecycleSplitter, DistillationColumn
from dme_reactor import HydraulicPump, ReversibleHeater, DMEEquilibriumReactor

class DigitalTwinEngine:
    """
    Motor do Xemelo Dixital de DME a partir de Biomasa.
    """
    def __init__(self):
        self.t = 0.0 # s
        self.dt = 2.0 # Paso de tempo base (s)
        
        # Inicializar todos os equipos
        # MÓDULO I
        self.feedstock = FeedstockManager()
        self.dryer = DryerRStoic()
        self.pyrolysis = PyrolysisRYield()
        self.gasifier = GasifierRGibbs()
        self.flash1 = FlashSeparator()
        
        # MÓDULO II
        self.compressor = MultiStageCompressor()
        self.r_meoh = MethanolPFR(N_nodes=5) # 5 nodos para velocidade en tempo real
        
        # Estado inicial do PFR (10 variables por nodo: 9 n_i, 1 T)
        self.state_meoh = np.zeros(self.r_meoh.N * 10)
        for k in range(self.r_meoh.N):
            # Condición inicial: arrincar con algo de N2 e algo de Metanol para acurtar o transitorio de arranque
            n_gas = (110e5 * self.r_meoh.V_cell) / (8.314 * 400.0) / 1000.0 # kmol
            self.state_meoh[k*10 + 3] = n_gas * 0.1 # 10% Metanol (Para que xa haxa caudal cara os separadores)
            self.state_meoh[k*10 + 6] = n_gas * 0.9 # 90% N2
            self.state_meoh[k*10 + 9] = 400.0 # K
        
        # MÓDULO III
        self.flash2_1 = DynamicFlash(T_set=50.0+273.15, P_set=10.0)
        self.splitter = RecycleSplitter(target_split=0.75)
        self.flash2_2 = DynamicFlash(T_set=40.0+273.15, P_set=2.6)
        # Columna S2-1 simplificada
        self.col_s2_1 = DistillationColumn(n_stages=10, feed_stage=5, RR=1.5, Q_reboiler=1e5, Q_condenser=1e5)
        
        # MÓDULO IV
        self.pump3 = HydraulicPump()
        self.heater3_1 = ReversibleHeater(target_T=100.0, name="HC-3.1")
        self.heater3_2 = ReversibleHeater(target_T=154.0, name="HC-3.2")
        self.r_dme = DMEEquilibriumReactor()
        
        # MÓDULO V
        self.col_d3_1 = DistillationColumn(n_stages=7, feed_stage=7, RR=2.0, Q_reboiler=23874.5*4.184, Q_condenser=35676.7*4.184)
        self.col_d3_2 = DistillationColumn(n_stages=15, feed_stage=11, RR=2.0, Q_reboiler=17383.8*4.184, Q_condenser=16457.0*4.184)
        
        # Correntes de reciclo con inventario
        # Orde: [H2, CO, CO2, CH3OH, H2O, DME, N2, O2, CH4]
        self.recycle_1_mol_h = np.zeros(9)
        self.recycle_2_mol_h = np.zeros(9)
        
        # Alimentación
        self.biomass_feed_kgh = 3000.0
        self.air_feed_kgh = 600.0
        
        # Inventario de presión do lazo de síntese
        self.synthesis_pressure = 10.0 # bar (para o reciclo)
        
        # Inicializar balance de masa por defecto
        self.last_mass_balance = {
            "in_biomass": self.biomass_feed_kgh,
            "in_air": self.air_feed_kgh,
            "out_dme": 0.0,
            "out_ash": 0.0,
            "out_water_f1": 0.0,
            "out_water_d32": 0.0,
            "out_purge_gas": 0.0,
        }

    def simulation_step(self):
        """
        Executa un paso temporal dt do xemelo dixital completo (Plantwide).
        """
        # =========================================================
        # MÓDULO I: Gasificación
        # =========================================================
        # 1. Secado
        biomass_dry, water_rel = self.dryer.compute(self.biomass_feed_kgh)
        
        # 2. Pirólise
        atoms_mol_h, ash_kgh, water_vapor = self.pyrolysis.compute(biomass_dry, water_rel, self.feedstock.get_composition())
        
        # 3. Gasificación
        # Dinámica de temperatura do gasificador
        self.gasifier.step_thermal(syngas_flow_mol_h=None, dt=self.dt)
        syngas_raw_mol_h = self.gasifier.compute_equilibrium(atoms_mol_h, self.air_feed_kgh)
        
        # 4. Separación Flash F1-1
        syngas_1_mol_h, liquid_1_mol_h, ash_out = self.flash1.compute(syngas_raw_mol_h, ash_kgh)
        self.flash1.step_dynamics(syngas_1_mol_h, liquid_1_mol_h, self.dt)
        
        # =========================================================
        # KPIs de Gasificación
        n_H2 = syngas_1_mol_h[0]
        n_CO = syngas_1_mol_h[1]
        n_CO2 = syngas_1_mol_h[2]
        
        syngas_M_ratio = 0.0
        if (n_CO + n_CO2) > 0:
            syngas_M_ratio = (n_H2 - n_CO2) / (n_CO + n_CO2)
            
        clean_syngas_kmol_h = np.sum(syngas_1_mol_h) / 1000.0

        # MÓDULO II: Compresión e Síntese de Metanol
        # =========================================================
        # Mesturar gas fresco co Reciclo I (Corrente R)
        syngas_mixed_mol_h = np.maximum(syngas_1_mol_h + self.recycle_1_mol_h, 0.0)
        
        # Calcular a presión real actual no primeiro nodo do reactor para o PID
        n_first = np.maximum(self.state_meoh[:9], 0.0)
        T_first = np.clip(self.state_meoh[9], 200.0, 1000.0)
        P_first_pa = (max(np.sum(n_first), 1e-10) * 1000.0 * 8.314 * T_first) / self.r_meoh.V_cell
        P_first_bar = P_first_pa / 1e5
        
        # Compresión
        # O PID de velocidade regula a presión
        rpm = self.compressor.step_dynamics(pv_P_out=P_first_bar, dt=self.dt)
        
        # Supoñendo que a presión se mantén preto do setpoint
        P_comp_out, T_comp_out, W_comp, n_stg, r_stg = self.compressor.compute(
            syngas_mixed_mol_h, P_in=self.flash1.P_current, T_in=40.0+273.15, target_P_out=self.compressor.P_out_setpoint)
            
        # PFR Metanol (Integración ríxida no tempo)
        t_span = (self.t, self.t + self.dt)
        feed_kmols_s = syngas_mixed_mol_h / 3600.0 / 1000.0
        
        def ode_func(t, y):
            return self.r_meoh.compute_derivatives(t, y, feed_kmols_s, T_feed=T_comp_out, u_valve=0.5, P_sink=10.0*1e5)
            
        sol = solve_ivp(ode_func, t_span, self.state_meoh, method='BDF')
        
        # Actualizar e acoutar o estado do reactor PFR para garantir estabilidade
        self.state_meoh = sol.y[:, -1]
        for k in range(self.r_meoh.N):
            self.state_meoh[k*10 : k*10 + 9] = np.maximum(self.state_meoh[k*10 : k*10 + 9], 0.0)
            self.state_meoh[k*10 + 9] = np.clip(self.state_meoh[k*10 + 9], 200.0, 1000.0)
        
        # Recuperar o fluxo de saída do último nodo do PFR (Corrente 9)
        # Moles no último nodo
        n_last = self.state_meoh[-10:-1]
        T_last = self.state_meoh[-1]
        sum_n_last = np.sum(n_last)
        if sum_n_last < 1e-10:
            y_last = np.zeros(9)
            y_last[6] = 1.0 # 100% N2
        else:
            y_last = n_last / sum_n_last
            
        P_last_pa = (np.sum(n_last) * 1000.0 * 8.314 * T_last) / self.r_meoh.V_cell
        P_last_bar = P_last_pa / 1e5
        
        # Fluxo saínte aproximado da válvula
        C_valve = 5e-7
        dP_valve = np.max([0.0, P_last_pa - 10.0*1e5])
        F_out_kmol_s = C_valve * 0.5 * np.sqrt(dP_valve)
        meoh_out_mol_h = np.maximum(y_last * F_out_kmol_s * 1000.0 * 3600.0, 0.0)
        
        # =========================================================
        # MÓDULO III: Tren de purificación de metanol e Reciclo I
        # =========================================================
        # V2-1 e Enfriador a 50C, 10 bar
        # Flash F2-1
        vapor_f21, liquid_f21 = self.flash2_1.compute_separation(meoh_out_mol_h)
        self.flash2_1.step_level(liquid_f21, self.dt)
        
        # Divisor de Reciclo SPT2-1 (Snowball mitigation via compositional purge)
        air_mol_h = (self.air_feed_kgh / 28.84) * 1000.0
        n2_in_mol_h = air_mol_h * 0.79
        
        recycle_gas, purge_gas, self.y_N2_current = self.splitter.split_dynamic_purge(vapor_f21, n2_in_mol_h)
        # Actualizar a variable de estado do reciclo
        self.recycle_1_mol_h = recycle_gas
        
        # Flash F2-2 (Líquido baixa de 10 a 2.6 bar)
        vapor_f22, liquid_f22 = self.flash2_2.compute_separation(liquid_f21)
        self.flash2_2.step_level(liquid_f22, self.dt)
        
        # Columna S2-1
        # Entra Corrente 15 (liquid_f22) + Corrente 27 (Reciclo II)
        feed_s21 = liquid_f22 + self.recycle_2_mol_h
        # Metanol purificado por colas (target 99.99%) - Simplificación: top e bot invertidos no modelo real?
        # A miúdo, lixeiros (gas) por cabeza, Metanol+Auga por colas cara D3-1. 
        # O purgado S2 por cabeza.
        top_s21, bot_s21 = self.col_s2_1.compute(feed_s21, purity_target_top=99.0, T_set_bot=150.0)
        
        # =========================================================
        # MÓDULO IV: Sección de Síntese de DME
        # =========================================================
        # P3-1
        P_dme_in, flow_dme_in = self.pump3.compute(bot_s21, P_in=2.6, P_out_set=15.1)
        
        # Quentadores reversibles H3-1 e H3-2
        T_h1, u_h1 = self.heater3_1.step_thermal(T_in=40.0+273.15, dt=self.dt)
        T_h2, u_h2 = self.heater3_2.step_thermal(T_in=T_h1, dt=self.dt)
        
        # Reactor R3-1 (Equilibrio)
        dme_out_mol_h = np.maximum(self.r_dme.compute_equilibrium(flow_dme_in), 0.0)
        
        # =========================================================
        # MÓDULO V: Purificación de DME e Control de Columnas
        # =========================================================
        # Columna D3-1 (Purificación de DME)
        # DME sae por cabeza
        top_d31, bot_d31 = self.col_d3_1.compute(dme_out_mol_h, purity_target_top=99.90, T_set_bot=155.0)
        
        # Columna D3-2 (Recuperación de Metanol e Reciclo II)
        # Entra bot_d31 (MeOH + H2O)
        # Por cabeza sae MeOH (Reciclo II), por colas H2O (residual)
        top_d32, bot_d32 = self.col_d3_2.compute(bot_d31, purity_target_top=97.32, T_set_bot=168.0)
        
        # Pechar o Reciclo II
        self.recycle_2_mol_h = top_d32
        
        # Gardar balance de masa (kg/h)
        # Pesos molares aprox: H2:2, CO:28, CO2:44, MeOH:32, H2O:18, DME:46, N2:28, O2:32, CH4:16
        MWs = np.array([2.016, 28.01, 44.01, 32.04, 18.015, 46.07, 28.013, 31.998, 16.04]) / 1000.0 # kg/mol
        
        DME_kgh = top_d31[5] * 46.07 / 1000.0
        purge_gas_kgh = np.sum(purge_gas * MWs)
        liquid_f1_kgh = np.sum(liquid_1_mol_h * MWs)
        bot_d32_kgh = np.sum(bot_d32 * MWs)
        ash_kgh_total = ash_out
        
        self.last_mass_balance = {
            "in_biomass": self.biomass_feed_kgh,
            "in_air": self.air_feed_kgh,
            "out_dme": DME_kgh,
            "out_ash": ash_kgh_total,
            "out_water_f1": liquid_f1_kgh,
            "out_water_d32": bot_d32_kgh,
            "out_purge_gas": purge_gas_kgh,
        }
        
        # Avanzar tempo
        self.t += self.dt
        
        # Resultados chave para a iteración
        results = {
            "time": self.t,
            "syngas_flow": np.sum(syngas_1_mol_h),
            "clean_syngas_kmol_h": clean_syngas_kmol_h,
            "syngas_M_ratio": syngas_M_ratio,
            "compressor_rpm": rpm,
            "compressor_stages": n_stg,
            "T_max_meoh": T_last - 273.15,
            "P_system": self.synthesis_pressure,
            "recycle_1_ratio": self.splitter.split_fraction,
            "y_N2": getattr(self, "y_N2_current", 0.0),
            "DME_production": top_d31[5] * 46.07 / 1000.0, # kg/h DME (MW aprox 46)
            "MeOH_recycle_2": np.sum(self.recycle_2_mol_h),
        }
        
        return results

if __name__ == "__main__":
    print("Iniciando Xemelo Dixital Plantwide...")
    engine = DigitalTwinEngine()
    
    # Simular 10 pasos para comprobar estabilidade
    for i in range(10):
        res = engine.simulation_step()
        print(f"t={res['time']:.1f}s | "
              f"Syngas={res['clean_syngas_kmol_h']:.1f} kmol/h (M={res['syngas_M_ratio']:.2f}) | "
              f"Comp={res['compressor_rpm']:.0f} RPM (stg={res['compressor_stages']}) | "
              f"T_Meoh={res['T_max_meoh']:.1f}C | "
              f"P_sys={res['P_system']:.2f} bar | "
              f"DME={res['DME_production']:.2f} kg/h")
