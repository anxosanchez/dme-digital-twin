# -*- coding: utf-8 -*-
"""
MÓDULO DE INTERACTION, PERTURBACIÓNS E CONTROL (control.py)
------------------------------------------------------------
Este módulo implementa os algoritmos de control industrial:
1. PID_Controller: Controlador PID discreto estándar con acción anti-windup (clamping).
2. DMC_Controller: Controlador Predictivo baseado en Modelo (MPC) tipo Dynamic Matrix Control (DMC).
"""

import numpy as np

class PIDController:
    """
    Controlador Proporcional-Integral-Derivativo (PID) discreto.
    Inclúe anti-windup por clamping (integración condicional) e límites de saída.
    """
    def __init__(self, Kp, Ki, Kd, Ts, u_min=0.0, u_max=1.0, action_type="reverse", name="PID"):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd
        self.Ts = Ts
        self.u_min = u_min
        self.u_max = u_max
        self.action_type = action_type # 'reverse' ou 'direct'
        self.name = name
        
        self.integral = 0.0
        self.prev_error = 0.0
        self.u = 0.0
        
    def reset(self, initial_u=0.0):
        self.integral = initial_u
        self.prev_error = 0.0
        self.u = initial_u

    def compute(self, sp, pv):
        """
        Calcula a acción de control u(k) a partir do setpoint (sp) e o process variable (pv).
        """
        if self.action_type == "direct":
            error = pv - sp
        else:
            error = sp - pv
            
        # Acción Proporcional
        P_term = self.Kp * error
        
        # Acción Derivativa (filtro simple para evitar amplificación do ruído)
        D_term = self.Kd * (error - self.prev_error) / self.Ts
        
        # Acción Integral provisional
        proposed_I = self.integral + self.Ki * error * self.Ts
        
        # Valor de saída sen saturar
        u_unsaturated = P_term + proposed_I + D_term
        
        # Aplicar límites á saída (Saturación)
        self.u = np.clip(u_unsaturated, self.u_min, self.u_max)
        
        # Anti-windup (Clamping): só integrar se o controlador non está saturado
        # ou se o erro e a saída teñen signos opostos (axuda a saír da saturación)
        is_saturated = (u_unsaturated != self.u)
        same_sign = np.sign(error) == np.sign(self.u - u_unsaturated)
        
        if not is_saturated or not same_sign:
            self.integral = proposed_I
            
        self.prev_error = error
        return self.u


class DMCController:
    """
    Controlador por Matriz Dinámica (DMC) - Modelo SISO (Single Input Single Output).
    Controla o reactor (CV = T_max) manipulando o refrixerante (MV = T_coolant).
    """
    def __init__(self, Kp=0.8, tau=150.0, theta=20.0, Ts=10.0, P=20, M=4, lambda_reg=1.0):
        """
        Kp: Ganancia do proceso (K/K)
        tau: Constante de tempo do proceso (s)
        theta: Tempo morto (s)
        Ts: Período de mostraxe (s)
        P: Horizonte de predición (número de intervalos)
        M: Horizonte de control (número de intervalos, M <= P)
        lambda_reg: Factor de penalización do movemento de control (regularización)
        """
        self.Kp = Kp
        self.tau = tau
        self.theta = theta
        self.Ts = Ts
        self.P = P
        self.M = M
        self.lambda_reg = lambda_reg
        
        # Xerar coeficientes da resposta ao chanzo s_i
        self.N_step = P + 20  # Número de coeficientes históricos
        self.s = np.zeros(self.N_step)
        for i in range(1, self.N_step + 1):
            t = i * Ts
            if t > theta:
                self.s[i-1] = Kp * (1.0 - np.exp(-(t - theta) / tau))
            else:
                self.s[i-1] = 0.0
                
        # Construír a Matriz Dinámica A (P x M)
        self.A = np.zeros((P, M))
        for j in range(M):
            self.A[j:, j] = self.s[:P-j]
            
        # Precalculo da matriz de ganancia de control: K_dmc = (A^T * A + lambda * I)^-1 * A^T
        ATA = np.dot(self.A.T, self.A)
        reg_matrix = lambda_reg * np.eye(M)
        self.inv_matrix = np.linalg.inv(ATA + reg_matrix)
        self.K_dmc = np.dot(self.inv_matrix, self.A.T)
        # Só aplicamos o primeiro movemento de control: K_first é a primeira fila de K_dmc
        self.K_first = self.K_dmc[0, :]
        
        # Inicialización de variables de estado do controlador
        self.past_dU = np.zeros(self.N_step) # Historial de cambios de control realizados
        self.u_prev = 0.0
        self.y_pred = np.zeros(self.P)       # Predición libre sen novos movementos

    def reset(self, initial_u):
        self.u_prev = initial_u
        self.past_dU = np.zeros(self.N_step)
        self.y_pred = np.zeros(self.P)

    def compute(self, sp, pv, u_min=450.0, u_max=600.0):
        """
        Calcula a nova acción de control u(k) dado o setpoint (sp) e a variable de proceso (pv).
        """
        # 1. Feedback e corrección do erro do modelo (Estimación do bias/disturbio d)
        # O bias é a diferenza entre o PV medido e a predición que fixemos no paso anterior para hoxe
        y_predicted_today = self.y_pred[0]
        bias = pv - y_predicted_today
        
        # 2. Actualizar a predición libre (free response) para o horizonte P
        # f_k = y_pred_k + bias
        # Onde desprazamos a predición un paso no tempo e actualizamos o final
        free_response = np.zeros(self.P)
        for i in range(self.P - 1):
            free_response[i] = self.y_pred[i+1] + bias
        # O último elemento estimase co valor asintótico
        free_response[self.P-1] = self.y_pred[self.P-1] + bias
        
        # 3. Calcular vector de erros respecto ao setpoint
        e_vector = sp - free_response
        
        # 4. Calcular o movemento de control óptimo (primeiro elemento do horizonte M)
        # dU = K_first * e_vector
        dU = np.dot(self.K_first, e_vector)
        
        # 5. Aplicar límites de velocidade e saturación absoluta ao MV
        # Máximo cambio por paso de tempo de 5 K
        dU_clipped = np.clip(dU, -5.0, 5.0)
        u_proposed = self.u_prev + dU_clipped
        u_final = np.clip(u_proposed, u_min, u_max)
        
        # Recalcular o dU realmente aplicado
        dU_applied = u_final - self.u_prev
        self.u_prev = u_final
        
        # 6. Actualizar o historial de dU
        # Desprazar cara atrás e insertar o novo dU_applied
        self.past_dU[1:] = self.past_dU[:-1]
        self.past_dU[0] = dU_applied
        
        # 7. Actualizar o vector de predicións para o próximo paso (k+1)
        # y_pred(k+1) = sum_{i} s_{i} * dU(k+1-i)
        # Calculamos as predicións para os seguintes P pasos
        for i in range(self.P):
            pred_val = 0.0
            for j in range(self.N_step):
                # O índice do dU debe ser relativo ao paso i-ésimo futuro
                # dU histórico está en past_dU: past_dU[0] = dU(k), past_dU[1] = dU(k-1)...
                idx_dU = j - i
                if idx_dU >= 0:
                    pred_val += self.s[j] * self.past_dU[idx_dU]
            self.y_pred[i] = pred_val
            
        return u_final


# Proba rápida de funcionamento do control
if __name__ == "__main__":
    # Test PID Reverse
    pid_rev = PIDController(Kp=0.5, Ki=0.01, Kd=0.1, Ts=1.0, action_type="reverse", u_min=0.0, u_max=100.0)
    print("TEST PID REVERSE:")
    print(f"PV=0, SP=100, u={pid_rev.compute(100, 0):.2f}")
    print(f"PV=50, SP=100, u={pid_rev.compute(100, 50):.2f}")
    
    # Test PID Direct
    pid_dir = PIDController(Kp=0.5, Ki=0.01, Kd=0.1, Ts=1.0, action_type="direct", u_min=0.0, u_max=100.0)
    print("\nTEST PID DIRECT:")
    print(f"PV=120, SP=100, u={pid_dir.compute(100, 120):.2f}")
    print(f"PV=50, SP=100, u={pid_dir.compute(100, 50):.2f}")
    
    # Test DMC
    dmc = DMCController(Kp=0.8, tau=100.0, theta=10.0, Ts=10.0, P=15, M=3, lambda_reg=0.5)
    dmc.reset(initial_u=513.15) # Inicializar a 240 C
    print("\nTEST DMC:")
    # Supoñemos que o proceso está en 523.15 C
    u_out = dmc.compute(sp=525.15, pv=523.15)
    print(f"PV=523.15, SP=525.15, u={u_out:.2f} K")
