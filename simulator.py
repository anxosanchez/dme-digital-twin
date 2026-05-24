# -*- coding: utf-8 -*-
"""
MOTOR CENTRAL E BUCLE DE SIMULACIÓN (simulator.py)
--------------------------------------------------
Este módulo implementa a clase principal `DMESimulator` que coordina:
- O estado dinámico do reactor de leito fixo.
- A integración numérica das ODEs ríxidas co solver BDF de SciPy.
- Os lazos de control regulando a presión (PID) e a temperatura (PID ou DMC).
- O rexistro histórico de variables de proceso (Series Temporais e perfís axiais).
"""

import numpy as np
from scipy.integrate import solve_ivp
from thermodynamics import R
from reactor import FixedBedReactor, CP_GAS, NU
from control import PIDController, DMCController

class DMESimulator:
    def __init__(self, N_nodes=10):
        # 1. Crear reactor e controladores
        self.reactor = FixedBedReactor(N_nodes=N_nodes, L=3.0, D_inner=0.5, mass_cat_total=800.0)
        self.N = N_nodes
        
        # Parámetros de simulación
        self.t = 0.0
        self.Ts = 2.0  # Paso de tempo do control e visualización (s)
        
        # Control de Presión (PID regulando a apertura de válvula)
        # sp = 50.0 bar, pv = P_outlet. Direct-acting (se P > 50, abrir válvula)
        self.pid_pressure = PIDController(
            Kp=-0.005, Ki=-0.002, Kd=-0.0002, Ts=self.Ts,
            u_min=0.05, u_max=1.0, name="PCV-101"
        )
        self.u_valve = 0.4  # Apertura inicial da válvula
        
        # Control de Temperatura (PID ou DMC regulando a temperatura do refrixerante T_coolant)
        # sp = T_setpoint, pv = T_max. Reverse-acting (se T > T_set, arrefriar)
        self.pid_temp = PIDController(
            Kp=1.2, Ki=0.008, Kd=0.05, Ts=self.Ts,
            u_min=450.0, u_max=570.0, name="TC-101 (PID)"
        )
        
        # DMC para temperatura: Horizonte predición P=25, control M=4, lambda=2.0
        self.dmc_temp = DMCController(
            Kp=0.75, tau=160.0, theta=16.0, Ts=self.Ts,
            P=25, M=4, lambda_reg=3.0
        )
        
        # Modo de control por defecto: 'PID' ou 'DMC' ou 'MANUAL'
        self.control_mode = "PID"
        self.T_setpoint = 533.15  # Setpoint de temperatura máxima: 260 °C
        self.T_coolant = 513.15   # Temperatura inicial do refrixerante: 240 °C (manual ou calculada por control)
        
        # Estado de alimentación e sink
        # Alimentación molar (mol/s): [H2, CO, CO2, CH3OH, H2O, DME]
        self.F_feed = np.array([6.0, 2.8, 0.4, 0.0, 0.0, 0.0]) # Total 9.2 mol/s
        self.T_feed = 523.15   # 250 C
        self.P_sink = 10.0 * 1e5  # 10 bar (presión downstream)
        
        # Perturbacións
        self.disturbance_flow_pct = 0.0  # Cambios na alimentación en %
        
        # Inicializar o estado do reactor e controladores
        self.reset_simulation()
        
    def reset_simulation(self):
        """
        Resetea a simulación ao estado inicial de 'Warmup'.
        """
        self.t = 0.0
        self.u_valve = 0.4
        self.T_coolant = 513.15
        
        # Estado inicial: gas con composición de entrada e temperatura uniforme de 250 C a 50 bar
        P_init = 50.0 * 1e5
        n_total_cell = (P_init * self.reactor.V_cell) / (R * self.T_feed)
        y_feed = self.F_feed / np.sum(self.F_feed)
        
        self.state = np.zeros(self.N * 7)
        for k in range(self.N):
            self.state[k*7 : k*7 + 6] = y_feed * n_total_cell
            self.state[k*7 + 6] = self.T_feed
            
        self.pid_pressure.reset(initial_u=self.u_valve)
        self.pid_temp.reset(initial_u=self.T_coolant)
        self.dmc_temp.reset(initial_u=self.T_coolant)
        
        # Historial de series temporais
        self.history = {
            "time": [],
            "T_max": [],
            "T_coolant": [],
            "P_outlet": [],
            "u_valve": [],
            "F_DME_outlet": [],
            "F_feed_total": [],
            "control_mode": []
        }
        
        # Historial de perfís espaciais (último paso de tempo)
        self.profile = {
            "z": np.linspace(0, self.reactor.L, self.N).tolist(),
            "T": [],
            "y_H2": [],
            "y_CO": [],
            "y_CO2": [],
            "y_MeOH": [],
            "y_H2O": [],
            "y_DME": [],
            "P": []
        }
        
        self.log_state()

    def log_state(self):
        """
        Rexistra o estado actual nas estruturas de datos de historial.
        """
        # Desempaquetar estado
        n = np.zeros((self.N, 6))
        T = np.zeros(self.N)
        for k in range(self.N):
            n[k, :] = self.state[k*7 : k*7 + 6]
            T[k] = self.state[k*7 + 6]
            
        n_total = np.sum(n, axis=1)
        y = n / n_total[:, np.newaxis]
        
        # Presión do último nodo (outlet) en bar
        P_est_last = (n_total[-1] * R * T[-1]) / self.reactor.V_cell
        Z_last, _ = self.reactor.eos.compute_mixture_properties(T[-1], P_est_last, y[-1, :])
        P_outlet_bar = (n_total[-1] * Z_last * R * T[-1]) / (self.reactor.V_cell * 1e5)
        
        # Moles e caudal molar de saída de DME
        C_valve = 8e-3
        P_last_pa = P_outlet_bar * 1e5
        dP_valve = max(0.0, P_last_pa - self.P_sink)
        F_total_valve = C_valve * self.u_valve * np.sqrt(dP_valve)
        F_DME_out = y[-1, 5] * F_total_valve
        
        # Rexistrar series temporais
        self.history["time"].append(self.t)
        self.history["T_max"].append(float(np.max(T) - 273.15)) # Celsius
        self.history["T_coolant"].append(float(self.T_coolant - 273.15)) # Celsius
        self.history["P_outlet"].append(float(P_outlet_bar))
        self.history["u_valve"].append(float(self.u_valve * 100)) # % de apertura
        self.history["F_DME_outlet"].append(float(F_DME_out))
        
        # Caudal de alimentación afectado por perturbacións
        F_actual_feed = self.F_feed * (1.0 + self.disturbance_flow_pct / 100.0)
        self.history["F_feed_total"].append(float(np.sum(F_actual_feed)))
        self.history["control_mode"].append(self.control_mode)
        
        # Limitar historial para aforrar memoria
        max_len = 300
        if len(self.history["time"]) > max_len:
            for key in self.history:
                self.history[key] = self.history[key][-max_len:]
                
        # Rexistrar perfile espacial actual
        self.profile["T"] = (T - 273.15).tolist() # Celsius
        self.profile["y_H2"] = y[:, 0].tolist()
        self.profile["y_CO"] = y[:, 1].tolist()
        self.profile["y_CO2"] = y[:, 2].tolist()
        self.profile["y_MeOH"] = y[:, 3].tolist()
        self.profile["y_H2O"] = y[:, 4].tolist()
        self.profile["y_DME"] = y[:, 5].tolist()
        
        # Presión por nodo en bar
        P_bar_nodes = []
        for k in range(self.N):
            P_est = (n_total[k] * R * T[k]) / self.reactor.V_cell
            Z_k, _ = self.reactor.eos.compute_mixture_properties(T[k], P_est, y[k, :])
            P_bar_nodes.append((n_total[k] * Z_k * R * T[k]) / (self.reactor.V_cell * 1e5))
        self.profile["P"] = P_bar_nodes

    def step(self):
        """
        Executa un paso de simulación temporal de tamaño Ts.
        - Lee sensores do reactor.
        - Executa algoritmos de control para calcular MV.
        - Integra as DAEs do reactor do instante t ao t+Ts.
        """
        # 1. Medir Variables de Proceso (PVs)
        # Obter perfís térmicos e de presión actuais
        n = np.zeros((self.N, 6))
        T = np.zeros(self.N)
        for k in range(self.N):
            n[k, :] = self.state[k*7 : k*7 + 6]
            T[k] = self.state[k*7 + 6]
            
        T_max = np.max(T)
        
        n_total = np.sum(n, axis=1)
        y = n / n_total[:, np.newaxis]
        P_est_last = (n_total[-1] * R * T[-1]) / self.reactor.V_cell
        Z_last, _ = self.reactor.eos.compute_mixture_properties(T[-1], P_est_last, y[-1, :])
        P_outlet_bar = (n_total[-1] * Z_last * R * T[-1]) / (self.reactor.V_cell * 1e5)
        
        # 2. Executar control de presión (PID)
        # Queremos P_outlet = 50.0 bar
        self.u_valve = self.pid_pressure.compute(sp=50.0, pv=P_outlet_bar)
        
        # 3. Executar control de temperatura (PID, DMC ou Manual)
        if self.control_mode == "PID":
            # sp = T_setpoint, pv = T_max
            # O PID de temperatura de control.py devolve a temperatura do coolant
            # Se T_max aumenta, queremos baixar T_coolant.
            # O PID de control.py foi configurado con erro = sp - pv.
            # sp - pv < 0 significa que pv > sp, e polo tanto a saída do PID debe diminuír.
            # Deseñamos a acción de control cun signo correspondente.
            self.T_coolant = self.pid_temp.compute(sp=self.T_setpoint, pv=T_max)
        elif self.control_mode == "DMC":
            # DMC calcula a temperatura de refrixeración
            self.T_coolant = self.dmc_temp.compute(sp=self.T_setpoint, pv=T_max, u_min=450.0, u_max=570.0)
        else:
            # Modo manual: non se modifica T_coolant fóra do control do usuario
            pass
            
        # 4. Aplicar perturbación de caudal se existe
        F_actual_feed = self.F_feed * (1.0 + self.disturbance_flow_pct / 100.0)
        
        # 5. Integración numérica ríxida co solver BDF de SciPy
        t_span = (self.t, self.t + self.Ts)
        
        # Ecuación diferencial encapsulada para o solver: dy/dt = f(t, y)
        def ode_func(t, y):
            return self.reactor.compute_derivatives(
                t, y,
                F_feed=F_actual_feed,
                T_feed=self.T_feed,
                T_coolant=self.T_coolant,
                u_valve=self.u_valve,
                P_sink=self.P_sink
            )
            
        # Chamada a solve_ivp
        sol = solve_ivp(
            ode_func,
            t_span,
            self.state,
            method='BDF', # Método axeitado para DAEs/ODEs ríxidas (stiff)
            rtol=1e-4,
            atol=1e-6
        )
        
        # Actualizar estado e tempo
        self.state = sol.y[:, -1]
        self.t += self.Ts
        
        # 6. Rexistrar datos
        self.log_state()


# Proba rápida da simulación e integración dinámicas
if __name__ == "__main__":
    sim = DMESimulator()
    print("Iniciando simulación de test...")
    print(f"Tempo inicial: {sim.t} s, T_max: {sim.history['T_max'][-1]:.2f} C, P_outlet: {sim.history['P_outlet'][-1]:.2f} bar")
    
    # Simular 5 pasos (10 segundos)
    for _ in range(5):
        sim.step()
        print(f"Tempo: {sim.t} s, T_max: {sim.history['T_max'][-1]:.2f} C, T_coolant: {sim.history['T_coolant'][-1]:.2f} C, P_outlet: {sim.history['P_outlet'][-1]:.2f} bar, Valve: {sim.history['u_valve'][-1]:.2f}%")
