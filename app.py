# -*- coding: utf-8 -*-
"""
MÓDULO VIII: ARQUITECTURA DO FRONTEND WEB (app.py)
--------------------------------------------------
Interface interactiva baseada en Streamlit tipo DCS (Distributed Control System).
"""

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import random
from engine import DigitalTwinEngine
from analytics import ProcessAnalytics

# ==========================================
# CONFIGURACIÓN DA PÁXINA E ESTILOS CSS
# ==========================================
st.set_page_config(page_title="DCS | DME Digital Twin", layout="wide", initial_sidebar_state="expanded")

CSS = """
<style>
    /* Tema escuro e DCS styling */
    .stApp {
        background-color: #0e1117;
        color: #c9d1d9;
    }
    
    /* Paneis con Glassmorphism */
    div[data-testid="stVerticalBlock"] > div[data-testid="stVerticalBlock"] {
        background: rgba(30, 34, 45, 0.4);
        border-radius: 8px;
        border: 1px solid rgba(255, 255, 255, 0.05);
        padding: 1.2rem;
    }
    
    /* Métricas destacadas */
    [data-testid="stMetricValue"] {
        font-size: 2.4rem !important;
        font-weight: 700 !important;
        color: #00ffcc !important;
    }
    [data-testid="stMetricLabel"] {
        font-size: 1.0rem !important;
        font-weight: 600 !important;
        color: #8b949e !important;
        text-transform: uppercase;
        letter-spacing: 1px;
    }
    
    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {
        gap: 24px;
        background-color: rgba(22, 27, 34, 0.8);
        padding: 0 15px;
        border-radius: 8px;
    }
    .stTabs [data-baseweb="tab"] {
        height: 50px;
        white-space: pre-wrap;
        background-color: transparent;
        border-radius: 4px 4px 0px 0px;
        gap: 1px;
        padding-top: 10px;
        padding-bottom: 10px;
        color: #c9d1d9;
    }
    
    /* Cabeceiras */
    h1, h2, h3 {
        color: #e6edf3 !important;
        font-weight: 600 !important;
    }
    
    /* Dataframes DCS */
    .dataframe {
        font-family: monospace;
        font-size: 0.95rem;
    }
</style>
"""
st.markdown(CSS, unsafe_allow_html=True)

# ==========================================
# CÁLCULO DINÁMICO DA TEMPERATURA (HUMIDADE)
# ==========================================
def calcular_resposta_humidade(tempo, modo_lazo, caudal_biomasa, queimadores_activos):
    # t é o vector de tempo de 0 a 60 min
    # Estado Inicial con Fallo (T cae de 750C a 650C)
    t = tempo
    if modo_lazo == "AUTOMÁTICO (AUTO)":
        # O PID automático satura e non pode recuperarse por si só do exceso de auga
        T = 750 - 100 * (1 - np.exp(-0.1 * t) * np.cos(0.05 * t))
    
    elif modo_lazo == "MANUAL (MAN)":
        # O cálculo depende do grao de corrección que faga o usuario co slider e o checkbox:
        # 1. Efecto do recorte de biomasa (reducir carga alivia o sumidoiro térmico)
        reduccion_carga = (3000 - caudal_biomasa) / 1200  # Normalizado de 0 a 1
        recuperacion_carga = 40 * reduccion_carga * (1 - np.exp(-0.12 * t))
        
        # 2. Efecto dos queimadores de emerxencia (+60C de achega térmica directa)
        if queimadores_activos:
            recuperacion_queimadores = 60 * (1 - np.exp(-0.2 * t))
        else:
            recuperacion_queimadores = 0
            
        # Ecuación final combinada
        T = 650 + recuperacion_carga + recuperacion_queimadores
        T = np.minimum(T, 750.0) # Teito en deseño
        
    return T

# ==========================================
# MOTOR MATEMÁTICO RADFRAC (D3-1)
# ==========================================
def calcular_mesh_torre(modo, R, F, DV):
    pratos = np.arange(1, 8)
    
    if modo == "Refluxo Total (D=0)":
        # Vapor interno segue o perfil nominal sen efecto do alimento F (F=0)
        V = 600 + 200 * np.exp(-0.25 * (pratos - 1))
        # Líquido L = V * (R / (R + 1)) * (1 + 0.1 * sin(prato))
        L = V * (R / (R + 1)) * (1 + 0.1 * np.sin(pratos))
        
        # Fraccións molares con separación máxima
        x_DME = 0.9999 * np.exp(-0.85 * (pratos - 1))
        x_MeOH = 0.8 * (1 - x_DME) * (pratos / 7)
        x_H2O = 1.0 - x_DME - x_MeOH
    else:
        # Operación Nominal (Produción)
        # O vapor inclúe o efeito do alimento F
        V = 600 + 200 * np.exp(-0.25 * (pratos - 1)) + (F * 0.4)
        # Líquido segue a ecuación acoplada
        L = V * (R / (R + 1)) * (1 + 0.1 * np.sin(pratos))
        
        # Fraccións molares nominais
        x_DME = 0.999 * np.exp(-0.5 * (pratos - 1))
        x_MeOH = 0.8 * (1 - x_DME) * (pratos / 7)
        x_H2O = 1.0 - x_DME - x_MeOH
        
    return pratos, L, V, x_DME, x_MeOH, x_H2O

# ==========================================
# INICIALIZACIÓN DO ESTADO E VARIABLES
# ==========================================
dt = 1.0  # Paso de integración global

if 'engine' not in st.session_state:
    st.session_state.engine = DigitalTwinEngine()
if 'analytics' not in st.session_state:
    st.session_state.analytics = ProcessAnalytics()
if 'is_running' not in st.session_state:
    st.session_state.is_running = False

# Variables de control do lazo TIC-01
if 'tic01_modo' not in st.session_state:
    st.session_state.tic01_modo = "AUTOMÁTICO (AUTO)"
if 'tic01_biomasa' not in st.session_state:
    st.session_state.tic01_biomasa = 3000
if 'tic01_queimadores' not in st.session_state:
    st.session_state.tic01_queimadores = False
if 'incident_start_time' not in st.session_state:
    st.session_state.incident_start_time = None

# Variables para a columna D3-1 (RadFrac)
if 'modo_operacion_col' not in st.session_state:
    st.session_state.modo_operacion_col = "Operación Nominal (Produción)"
if 'plates_F_new' not in st.session_state:
    st.session_state.plates_F_new = 689
if 'plates_R_nominal' not in st.session_state:
    st.session_state.plates_R_nominal = 2.0
if 'plates_DV_nominal' not in st.session_state:
    st.session_state.plates_DV_nominal = 0.33

# Variables protexidas para métricas e UI
if 'caudal_dme' not in st.session_state:
    st.session_state.caudal_dme = 0.00
if 'hotspot_t' not in st.session_state:
    st.session_state.hotspot_t = 150.0
if 'presion_sistema' not in st.session_state:
    st.session_state.presion_sistema = 1.0
if 'nitroxeno_acumulado' not in st.session_state:
    st.session_state.nitroxeno_acumulado = 0.0

# Variables de simulación P&ID
if 'pid_sim' not in st.session_state:
    st.session_state.pid_sim = {
        'time': [],
        'T_gasifier': [],
        'R_recycle': []
    }

# ==========================================
# CÁLCULO DINÁMICO DE ACOPLAMENTO RADFRAC (D3-1)
# ==========================================
modo_dest = st.session_state.modo_operacion_col
F_dest = st.session_state.plates_F_new

if modo_dest == "Operación Nominal (Produción)":
    R_dest = st.session_state.plates_R_nominal
    DV_dest = st.session_state.plates_DV_nominal
else:
    R_dest = 9999.0
    DV_dest = 0.0

# Executar a emulación MESH da columna
pratos_mesh, L_mesh, V_mesh, x_DME_mesh, x_MeOH_mesh, x_H2O_mesh = calcular_mesh_torre(modo_dest, R_dest, F_dest, DV_dest)

# Condición de corrente de saída condicional
modo_refluxo = modo_dest
D_V_split = DV_dest

if modo_refluxo == "Refluxo Total (D=0)":
    caudal_DME_produto = 0.0
    pureza_DME_comercial = 0.0
else:
    caudal_DME_produto = V_mesh[0] * D_V_split
    pureza_DME_comercial = x_DME_mesh[0] * 100.0

if not st.session_state.is_running:
    caudal_DME_produto = 0.0

engine = st.session_state.engine
analytics = st.session_state.analytics

# ==========================================
# BARRA LATERAL (SIDEBAR)
# ==========================================
st.sidebar.title("🎛️ DCS Control")
st.sidebar.markdown("---")

# Alerta global de Refluxo Total na barra lateral
if st.session_state.modo_operacion_col == "Refluxo Total (D=0)":
    st.sidebar.error("🚨 **PRODUCIÓN DETIDA:** A columna D3-1 está en **Refluxo Total**. A planta atópase en ciclo pechado sen saída comercial.")

btn_run = st.sidebar.button("▶️ Executar Simulación", use_container_width=True)
if btn_run:
    st.session_state.is_running = not st.session_state.is_running

# Manter a variable btn_run conectada ao estado persistente para o resto do script
btn_run = st.session_state.is_running

# Inicializar incidente no session_state
if 'active_incident' not in st.session_state:
    st.session_state.active_incident = "Operación Nominal (Segura)"

# Controis SCADA dinámicos para o lazo TIC-01 se o escenario é Biomasa Húmida
if st.session_state.active_incident == "Inxección de Biomasa Húmida (Chuvia no Silo)":
    st.sidebar.markdown("### 🎛️ Panel Lazo TIC-01 (Biomasa Húmida)")
    
    st.session_state.tic01_modo = st.sidebar.radio(
        "Modo do Lazo TIC-01:",
        ["AUTOMÁTICO (AUTO)", "MANUAL (MAN)"],
        index=0 if st.session_state.tic01_modo == "AUTOMÁTICO (AUTO)" else 1
    )
    
    is_manual = (st.session_state.tic01_modo == "MANUAL (MAN)")
    
    st.session_state.tic01_biomasa = st.sidebar.slider(
        "Alimentación de Biomasa (kg/h):",
        min_value=1500,
        max_value=3000,
        value=st.session_state.tic01_biomasa,
        step=50,
        disabled=not is_manual
    )
    
    st.session_state.tic01_queimadores = st.sidebar.toggle(
        "Activar Queimadores Auxiliares de Gasóleo (Emerxencia)",
        value=st.session_state.tic01_queimadores
    )
    st.sidebar.markdown("---")

# Lazo de asignación de variables de planta baseados na simulación (btn_run e active_incident)
if btn_run:
    if st.session_state.active_incident == "Operación Nominal (Segura)":
        st.session_state.caudal_dme = 438.47
        st.session_state.hotspot_t = 267.0
        st.session_state.presion_sistema = 110.0
    elif st.session_state.active_incident == "Inxección de Biomasa Húmida (Chuvia no Silo)":
        # Cálculo dinámico basado en temperatura
        if st.session_state.incident_start_time is None:
            t_gas_base = 650.0
        else:
            elapsed_time_min = (engine.t - st.session_state.incident_start_time) / 60.0
            t_gas_base = calcular_resposta_humidade(
                elapsed_time_min,
                st.session_state.tic01_modo,
                st.session_state.tic01_biomasa,
                st.session_state.tic01_queimadores
            )
        fraction = np.clip((t_gas_base - 650.0) / 100.0, 0.0, 1.0)
        st.session_state.caudal_dme = 241.16 + (438.47 - 241.16) * fraction
        st.session_state.hotspot_t = 210.5 + (267.0 - 210.5) * fraction
        st.session_state.presion_sistema = 110.0
    elif st.session_state.active_incident == "Bloqueo de Válvula de Purga (Efecto Bóla de Neve de N2)":
        st.session_state.caudal_dme = 21.92
        st.session_state.hotspot_t = 162.1
        st.session_state.presion_sistema = 110.0
    elif st.session_state.active_incident == "Inundación por Pico de Caudal na Torre D3-1":
        st.session_state.caudal_dme = 412.16
        st.session_state.hotspot_t = 267.0
        st.session_state.presion_sistema = 110.0
else:
    st.session_state.caudal_dme = 0.00
    st.session_state.hotspot_t = 150.0
    st.session_state.presion_sistema = 1.0

if st.sidebar.button("🔄 RESET", use_container_width=True):
    st.session_state.engine = DigitalTwinEngine()
    st.session_state.is_running = False
    st.session_state.active_incident = "Operación Nominal (Segura)"
    st.session_state.incident_start_time = None
    st.session_state.tic01_modo = "AUTOMÁTICO (AUTO)"
    st.session_state.tic01_biomasa = 3000
    st.session_state.tic01_queimadores = False
    st.session_state.caudal_dme = 0.00
    st.session_state.hotspot_t = 150.0
    st.session_state.presion_sistema = 1.0
    st.session_state.nitroxeno_acumulado = 0.0
    st.session_state.pid_sim = {'time': [], 'T_gasifier': [], 'R_recycle': []}
    st.rerun()

st.sidebar.markdown("---")
st.sidebar.subheader("🌲 Feedstock Matrix (%)")

# Porcentaxes reais de piñeiro como base
c_val = st.sidebar.slider("Carbono (C)", 0.0, 100.0, 43.45, 0.01)
h_val = st.sidebar.slider("Hidróxeno (H)", 0.0, 100.0, 4.24, 0.01)
o_val = st.sidebar.slider("Osíxeno (O)", 0.0, 100.0, 39.10, 0.01)
n_val = st.sidebar.slider("Nitróxeno (N)", 0.0, 100.0, 0.72, 0.01)
s_val = st.sidebar.slider("Xofre (S)", 0.0, 100.0, 0.40, 0.01)
ash_val = st.sidebar.slider("Cinzas", 0.0, 100.0, 12.09, 0.01)

sum_comp = c_val + h_val + o_val + n_val + s_val + ash_val

# Normalización dinámica
if sum_comp > 0:
    c_norm = (c_val / sum_comp) * 100.0
    h_norm = (h_val / sum_comp) * 100.0
    o_norm = (o_val / sum_comp) * 100.0
    n_norm = (n_val / sum_comp) * 100.0
    s_norm = (s_val / sum_comp) * 100.0
    ash_norm = (ash_val / sum_comp) * 100.0
else:
    c_norm, h_norm, o_norm, n_norm, s_norm, ash_norm = 100.0, 0.0, 0.0, 0.0, 0.0, 0.0

if abs(sum_comp - 100.0) > 0.01:
    st.sidebar.warning(f"Suma orixinal: {sum_comp:.1f}%. A matriz foi auto-normalizada ao 100.0%.")
    st.sidebar.info(f"Matriz real: C:{c_norm:.1f}% | H:{h_norm:.1f}% | O:{o_norm:.1f}% | N:{n_norm:.2f}% | S:{s_norm:.2f}% | Cinzas:{ash_norm:.1f}%")

if st.sidebar.button("💾 Aplicar Feedstock", use_container_width=True):
    engine.feedstock.set_custom_composition(c_norm/100.0, h_norm/100.0, o_norm/100.0, n_norm/100.0, s_norm/100.0, ash_norm/100.0)
    st.sidebar.success("Composición normalizada rexistrada no solver.")

# Lóxica dinámica de simulación
if btn_run:
    engine.active_incident = st.session_state.active_incident
    
    t_gas_noise = random.uniform(-1.5, 1.5)
    
    if st.session_state.active_incident == "Inxección de Biomasa Húmida (Chuvia no Silo)":
        if st.session_state.incident_start_time is None:
            st.session_state.incident_start_time = engine.t
        elapsed_time_min = (engine.t - st.session_state.incident_start_time) / 60.0
        
        t_gas_base = calcular_resposta_humidade(
            elapsed_time_min,
            st.session_state.tic01_modo,
            st.session_state.tic01_biomasa,
            st.session_state.tic01_queimadores
        )
        t_gas = t_gas_base + t_gas_noise
        
        # Set the override so engine respects it
        engine.gasifier_T_override = t_gas_base + 273.15
        
        if st.session_state.tic01_modo == "MANUAL (MAN)":
            engine.biomass_feed_kgh = float(st.session_state.tic01_biomasa)
        else:
            engine.biomass_feed_kgh = 3000.0
            
        fraction = np.clip((t_gas_base - 650.0) / 100.0, 0.0, 1.0)
        st.session_state.caudal_dme = 241.16 + (438.47 - 241.16) * fraction
        st.session_state.hotspot_t = 210.5 + (267.0 - 210.5) * fraction
        st.session_state.presion_sistema = 110.0
    else:
        st.session_state.incident_start_time = None
        engine.gasifier_T_override = None
        t_gas_base = 750.0
        t_gas = t_gas_base + t_gas_noise
        engine.biomass_feed_kgh = 3000.0

    res = engine.simulation_step()
    st.session_state.nitroxeno_acumulado = res.get("y_N2", 0.0) * 100.0
    
    # Alimentar o histórico do P&ID
    current_time = engine.t
    
    # Setpoint FIC-02 = 75.0%
    fic_noise = random.uniform(-0.5, 0.5)
    fic_val = 75.0 + fic_noise
    
    st.session_state.pid_sim['time'].append(current_time)
    st.session_state.pid_sim['T_gasifier'].append(t_gas)
    st.session_state.pid_sim['R_recycle'].append(fic_val)
    
    # Manter só os últimos 50 puntos
    if len(st.session_state.pid_sim['time']) > 50:
        st.session_state.pid_sim['time'] = st.session_state.pid_sim['time'][-50:]
        st.session_state.pid_sim['T_gasifier'] = st.session_state.pid_sim['T_gasifier'][-50:]
        st.session_state.pid_sim['R_recycle'] = st.session_state.pid_sim['R_recycle'][-50:]

    pass # st.rerun() moved to the bottom of the script


# ==========================================
# CABECEIRA E KPIs
# ==========================================
st.title("🏭 Xemelo Dixital de Síntese de DME (DCS)")

# Inicialización segura para o primeiro render antes de usar os botóns
if 'caudal_dme' not in st.session_state:
    st.session_state.caudal_dme = 0.00
if 'hotspot_t' not in st.session_state:
    st.session_state.hotspot_t = 150.0
if 'presion_sistema' not in st.session_state:
    st.session_state.presion_sistema = 1.0
if 'nitroxeno_acumulado' not in st.session_state:
    st.session_state.nitroxeno_acumulado = 0.0

kpi_dme_delta = "Activo" if btn_run else "Parado"

col1, col2, col3, col4 = st.columns(4)
col1.metric("Caudal de DME", f"{caudal_DME_produto:.2f} kg/h", delta=kpi_dme_delta)
col2.metric("Hotspot Metanol", f"{st.session_state.hotspot_t:.1f} °C", delta=f"{st.session_state.hotspot_t - 150.0:.1f} °C" if st.session_state.hotspot_t > 150 else None)
col3.metric("K2-1 Descarga", f"{st.session_state.presion_sistema:.1f} bar", delta=f"{st.session_state.presion_sistema - 1.0:.1f} bar" if st.session_state.presion_sistema > 1 else None)
col4.metric("Acumulación N2", f"{st.session_state.nitroxeno_acumulado:.2f} %", delta="SP: 18.54 %" if btn_run else None)

st.markdown("---")

# ==========================================
# PESTANAS DE VISUALIZACIÓN
# ==========================================
tab_scada, tab_pid, tab_analytics, tab_incidents, tab_plates = st.tabs([
    "📊 SCADA Overview", 
    "🎛️ Control P&ID (Vivo)", 
    "📈 Motor Analítico e Balances",
    "🚨 Simulación de Incidencias",
    "📊 Perfil de Pratos (D3-1)"
])

# ----------------- TAB 1: SCADA -----------------
with tab_scada:
    # Variables dinámicas de visualización SCADA baseadas no incidente e no estado de marcha
    t_gasifier_scada = 750.0
    caudal_syngas_scada = 3057.81
    purity_dme_scada = pureza_DME_comercial
    delta_p_d31_scada = 20.0 # mbar nominal
    
    if btn_run:
        if st.session_state.active_incident == "Inxección de Biomasa Húmida (Chuvia no Silo)":
            if len(st.session_state.pid_sim['T_gasifier']) > 0:
                t_gasifier_scada = st.session_state.pid_sim['T_gasifier'][-1]
            else:
                t_gasifier_scada = 650.0
            fraction = np.clip((t_gasifier_scada - 650.0) / 100.0, 0.0, 1.0)
            caudal_syngas_scada = 1987.58 + (3057.81 - 1987.58) * fraction
        elif st.session_state.active_incident == "Inundación por Pico de Caudal na Torre D3-1":
            purity_dme_scada = 94.20
            delta_p_d31_scada = 80.0 # 400% aumento
            
    col_left, col_right = st.columns(2)
    
    with col_left:
        st.markdown("### 🔥 Área Quente (Atm)")
        st.markdown(f"""
        **Sección de Gasificación**
        - Temperatura Reactor (RGibbs): `{t_gasifier_scada:.1f} °C`
        - Presión de operación: `1.0 bar`
        - Caudal de Syngas Limpo (F1-1): `{caudal_syngas_scada:.2f} kg/h`
        
        **Sección de Compresión (K2-1)**
        - Presión de Aspiración: `1.0 bar`
        - Presión de Descarga: `{st.session_state.presion_sistema:.1f} bar`
        - Etapas Activas: `5` (Ratio ~2.55)
        - Velocidade (RPM): `{'2500' if btn_run else '0'}`
        """)
        
    with col_right:
        st.markdown("### 🧪 Área de Alta Presión")
        st.markdown(f"""
        **Sección de Metanol (R2-1)**
        - Presión de Síntese PFR: `110.0 bar`
        - Lazo de Reciclo I (FIC-02): `75.0 %`
        
        **Sección de DME (R3-1 / D3-1)**
        - Temp. Entrada REquil (H3-2): `154.0 °C`
        - Presión Reactor DME: `14.7 bar`
        - Pureza DME Cabeza D3-1: `{purity_dme_scada:.2f} %`
        - Caída de Presión en D3-1 (Delta P): `{delta_p_d31_scada:.1f} mbar`
        """)

# ----------------- TAB 2: P&ID -----------------
with tab_pid:
    st.markdown("### Lazos de Control Dinámico (Anti-Windup)")
    
    if btn_run and len(st.session_state.pid_sim['time']) > 0:
        col_chart1, col_chart2 = st.columns(2)
        
        with col_chart1:
            st.markdown("**TIC-01: Temperatura Gasificador** (SP: 750 °C)")
            df_t = pd.DataFrame({
                'Tempo (s)': st.session_state.pid_sim['time'],
                'Temperatura (°C)': st.session_state.pid_sim['T_gasifier']
            }).set_index('Tempo (s)')
            st.line_chart(df_t, color="#ff4b4b")
            
        with col_chart2:
            st.markdown("**FIC-02: Válvula Reciclo I** (SP: 75.0 %)")
            df_r = pd.DataFrame({
                'Tempo (s)': st.session_state.pid_sim['time'],
                'Reciclo (%)': st.session_state.pid_sim['R_recycle']
            }).set_index('Tempo (s)')
            st.line_chart(df_r, color="#00ffcc")
    else:
        st.info("Inicia a Planta dende o DCS Control (Barra Lateral) para monitorizar os lazos PID en vivo.")

# ----------------- TAB 3: ANALÍTICA -----------------
with tab_analytics:
    st.markdown("### Informes Avanzados de Enxeñaría (En Vivo)")
    
    # 1. Gráficas lado a lado
    col_chart_left, col_chart_right = st.columns(2)
    
    with col_chart_left:
        st.markdown("#### 1. Sensibilidade Térmica do Gasificador (RGibbs)")
        
        comp = engine.feedstock.get_composition()
        comp_hash = hash(frozenset(comp.items()))
        
        if 'gasifier_sensitivity_fig' not in st.session_state or st.session_state.get('last_comp_hash') != comp_hash:
            if 'gasifier_sensitivity_fig' in st.session_state:
                plt.close(st.session_state.gasifier_sensitivity_fig)
            fig_sens = analytics.plot_gasifier_sensitivity(engine)
            st.session_state.gasifier_sensitivity_fig = fig_sens
            st.session_state.last_comp_hash = comp_hash
            
        st.pyplot(st.session_state.gasifier_sensitivity_fig)
        
    with col_chart_right:
        st.markdown("#### 2. Perfil Axial no Reactor PFR de Metanol (R2-1)")
        fig_pfr = analytics.plot_pfr_profiles(engine)
        st.pyplot(fig_pfr)
        plt.close(fig_pfr)
        
    st.markdown("---")
    
    # 2. Táboas de balances de materia e elemental lado a lado
    st.markdown("### 📊 Balances Operativos do Proceso (Dinámicos)")
    
    col_table_left, col_table_right = st.columns(2)
    
    # Obter balances dinámicos
    df_mass, df_elem = analytics.generate_mass_balance_report(engine)
    
    with col_table_left:
        st.markdown("#### Balance de Materia Global por Correntes")
        st.table(df_mass.set_index("Corrente"))
        
        total_in = df_mass[df_mass["Tipo"] == "Entrada"]["Caudal (kg/h)"].sum()
        total_out = df_mass[df_mass["Tipo"] == "Saída"]["Caudal (kg/h)"].sum() + df_mass[df_mass["Tipo"] == "Inventario"]["Caudal (kg/h)"].sum()
        desviacion = abs(total_in - total_out)
        desv_pct = (desviacion / total_in) * 100.0 if total_in > 0 else 0.0
        st.success(f"✔️ **Balance pechado.** Desviación: {desviacion:.4f} kg/h ({desv_pct:.2f}%)")
        
    with col_table_right:
        st.markdown("#### Balance Elemental por Especies")
        st.table(df_elem.set_index("Especie Elemental"))

# ----------------- TAB 4: INCIDENCIAS -----------------
with tab_incidents:
    st.markdown("### 🚨 Simulación de Incidencias e Perturbacións Críticas")
    st.markdown("Use esta sección para estresar os lazos de control da planta e ver a resposta dinámica en tempo real.")
    
    # Selector de incidencias
    opcion_incidencia = st.selectbox(
        "Seleccione a perturbación ou escenario a simular:",
        [
            "Operación Nominal (Segura)",
            "Inxección de Biomasa Húmida (Chuvia no Silo)",
            "Bloqueo de Válvula de Purga (Efecto Bóla de Neve de N2)",
            "Inundación por Pico de Caudal na Torre D3-1"
        ],
        index=[
            "Operación Nominal (Segura)",
            "Inxección de Biomasa Húmida (Chuvia no Silo)",
            "Bloqueo de Válvula de Purga (Efecto Bóla de Neve de N2)",
            "Inundación por Pico de Caudal na Torre D3-1"
        ].index(st.session_state.active_incident)
    )
    
    if opcion_incidencia != st.session_state.active_incident:
        st.session_state.active_incident = opcion_incidencia
        st.rerun()
        
    st.markdown("---")
    
    # Bloques visuais de feedback e accións
    if st.session_state.active_incident == "Operación Nominal (Segura)":
        st.success("🟢 **Operación Nominal (Segura)**: Planta operando en rango seguro. Todos os sistemas de control en automático (AUTO).")
        st.info("ℹ️ **Lazos PID activos**: Controladores sintonizados en modo AUTO mantendo setpoints nominais sen oscilacións significativas.")
        st.markdown("""
        <div style="background-color: #1e293b; border-left: 5px solid #10b981; padding: 15px; border-radius: 4px;">
            <strong style="color: #10b981;">✅ Acción Humana Requirida:</strong><br>
            Sen accións de emerxencia requiridas. Monitorizar variables a través do SCADA e realizar controis analíticos rutineiros cada 2 horas.
        </div>
        """, unsafe_allow_html=True)
        
    elif st.session_state.active_incident == "Inxección de Biomasa Húmida (Chuvia no Silo)":
        st.warning("⚠️ **Síntomas Físicos na Planta**: Caída drástica da temperatura do gasificador por evaporación endotérmica da humidade. Redución do caudal de syngas limpo e desprazamento da relación molar M. Risco de formação de alcatráns no leito do gasificador.")
        st.info("ℹ️ **Acción PID Automática**: O lazo TIC-01 (temperatura do gasificador) detecta a caída térmica e abre ao 100% a válvula de aire secundario/combustión, pero a capacidade de recuperación térmica está saturada pola alta humidade.")
        
        # Alertas de seguridade dinámicas baseadas no protocolo
        is_safe = (
            st.session_state.tic01_modo == "MANUAL (MAN)" and 
            st.session_state.tic01_biomasa < 2000 and 
            st.session_state.tic01_queimadores
        )
        
        if is_safe:
            st.success("🟢 SEGURO: Perfil térmico recuperado por enriba da zona de alcatráns. Planta estabilizada.")
        else:
            st.error("🚨 CRÍTICO: Temperatura inferior a 700°C. Iníciase a xeración de alcatráns. Risco de taponamento do leito de biomasa!")
            
        st.markdown("""
        <div style="background-color: #1e293b; border-left: 5px solid #f59e0b; padding: 15px; border-radius: 4px;">
            <strong style="color: #f59e0b;">🛠️ Acción Humana Requirida:</strong><br>
            Pasar o lazo TIC-01 a modo Manual (MAN). Reducir a alimentación de biomasa de 3000 kg/h a 1800 kg/h. Activar os queimadores auxiliares de gasóleo de emerxencia para recuperar o perfil térmico superior a 700°C e evitar a xeración de alcatráns que poidan taponar o leito.
        </div>
        """, unsafe_allow_html=True)
        
    elif st.session_state.active_incident == "Bloqueo de Válvula de Purga (Efecto Bóla de Neve de N2)":
        st.error("🚨 **Síntomas Físicos na Planta**: Acumulación descontrolada de Nitróxeno (inerte) no lazo de reciclo I. A fração molar de N2 elévase ata o 42%, asfixiando as presións parciais de H2 e CO. Paradas parciais de reacción no PFR de Metanol por asfixia termodinámica. Perda total de rendemento.")
        st.info("ℹ️ **Acción PID Automática**: O controlador de purga tenta compensar abrindo a válvula de control de purga ao límite físico (100% OP), pero o bloqueo mecánico augas abaixo impide a saída de gas inerte, inutilizando o control de snowball.")
        st.markdown("""
        <div style="background-color: #1e293b; border-left: 5px solid #ef4444; padding: 15px; border-radius: 4px;">
            <strong style="color: #ef4444;">🛠️ Acción Humana Requirida:</strong><br>
            Activar inmediatamente a liña de Bypass manual da purga cara ao facho (flare). Se a fracción de N2 supera o 40% durante máis de 5 minutos, proceder a Parada de Emerxencia (ESD-1) do lazo de síntese para protexer o catalizador de metanol contra a desactivación térmica.
        </div>
        """, unsafe_allow_html=True)
        
    elif st.session_state.active_incident == "Inundación por Pico de Caudal na Torre D3-1":
        st.error("🚨 **Síntomas Físicos na Planta**: Aumento brusco da caída de presión (Delta P) na columna D3-1 ata os 80 mbar (400% do nominal). Perda do perfil de temperaturas, inundación física de pratos e arrastre de fração pesada (auga e metanol) por cabeza de columna, derrubando a pureza de DME do 99.9% ao 94.2%.")
        st.info("ℹ️ **Acción PID Automática**: O lazo TIC-03 detecta alteración de temperatura no prato de control e incrementa ao máximo o caudal de refluxo (L) cara á cabeça da torre, o que agrava a inundación debido ao pico de caudal volumétrico de líquido.")
        st.markdown("""
        <div style="background-color: #1e293b; border-left: 5px solid #ef4444; padding: 15px; border-radius: 4px;">
            <strong style="color: #ef4444;">🛠️ Acción Humana Requirida:</strong><br>
            Reducir de inmediato o refluxo da torre. Diminuír a potencia de calefacción no calderín (Q_reboiler) para baixar a velocidade do vapor ascensional. Se a inundación persiste, pasar a torre a modo de refluxo total para estabilizar os pratos e evitar o arrastre continuo.
        </div>
        """, unsafe_allow_html=True)

# ----------------- TAB 5: PERFIL DE PRATOS (D3-1) -----------------
with tab_plates:
    st.markdown("### 📊 Perfil de Pratos Teóricos en D3-1 (RadFrac Emulation)")
    st.markdown("Esta sección modela o comportamento hidrodinámico e as fraccións molares líquidas prato a prato en tempo real.")
    
    # 1. Controis interactivos (sliders) en columnas superiores
    col_mode, col_feed = st.columns(2)
    with col_mode:
        modo_operacion = st.radio(
            "Modo de Operación da Columna:",
            ["Operación Nominal (Produción)", "Refluxo Total (D=0)"],
            horizontal=True
        )
    with col_feed:
        F_val = st.slider(
            "Caudal de Alimentación dende R3-1 (kg/h)",
            min_value=300,
            max_value=1000,
            value=689,
            step=10,
            key="plates_F_new"
        )
        
    st.markdown("---")
    
    # Colocar os sliders de fraccionamiento condicionales
    col_r, col_dv = st.columns(2)
    
    if modo_operacion == "Operación Nominal (Produción)":
        with col_r:
            R_val = st.slider(
                "Relación de Refluxo Exterior (R = L/D)",
                min_value=0.5,
                max_value=5.0,
                value=2.0,
                step=0.1,
                key="plates_R_nominal"
            )
        with col_dv:
            DV_val = st.slider(
                "Relación Destilado / Vapor de Cabezas (D/V)",
                min_value=0.1,
                max_value=0.9,
                value=0.33,
                step=0.01,
                key="plates_DV_nominal"
            )
    else:
        # Refluxo Total
        R_val = 9999.0
        DV_val = 0.0
        with col_r:
            st.info("ℹ️ **Refluxo Total Activo:** R fixado en ∞ (9999.0)")
        with col_dv:
            st.info("ℹ️ **Destilado Pechado:** D/V fixado en 0.0")

    # 2. Motor de cálculo MESH integrado
    pratos, L, V, x_DME, x_MeOH, x_H2O = calcular_mesh_torre(modo_operacion, R_val, F_val, DV_val)
    
    # Crear nomes dos pratos para o eixe y da gráfica
    prato_names = [f"Prato {i}" if i not in [1, 7] else (f"Prato 1 (Condensador)" if i == 1 else f"Prato 7 (Refervedor)") for i in pratos]
    
    # 3. Requisitos gráficos (visualización SCADA)
    col_g1, col_g2 = st.columns(2)
    
    with col_g1:
        st.markdown("#### 🌊 Hidrodinámica Interna (Líquido vs Vapor)")
        df_hydro = pd.DataFrame({
            "Prato": prato_names,
            "Líquido (L)": L,
            "Vapor (V)": V
        }).set_index("Prato")
        # Invertir para que o prato 1 quede arriba e o 7 abaixo
        df_hydro_rev = df_hydro.iloc[::-1]
        st.bar_chart(df_hydro_rev, horizontal=True, stack=False)
        
    with col_g2:
        st.markdown("#### 🧪 Perfil de Fraccións Molares (x_i)")
        df_chem = pd.DataFrame({
            "Prato": prato_names,
            "DME (x_DME)": x_DME,
            "Metanol (x_MeOH)": x_MeOH,
            "Auga (x_H2O)": x_H2O
        }).set_index("Prato")
        # Invertir para que o prato 1 quede arriba e o 7 abaixo
        df_chem_rev = df_chem.iloc[::-1]
        st.bar_chart(df_chem_rev, horizontal=True, stack=True)
        
    # 4. Matriz de datos en tempo real
    st.markdown("#### 📋 Matriz MESH de Datos Analíticos")
    df_mesh_table = pd.DataFrame({
        "Etapa / Prato": prato_names,
        "Caudal Líquido L (kg/h)": L,
        "Caudal Vapor V (kg/h)": V,
        "Frac. Molar DME (x_DME)": x_DME,
        "Frac. Molar MeOH (x_MeOH)": x_MeOH,
        "Frac. Molar Auga (x_H2O)": x_H2O
    }).set_index("Etapa / Prato")
    
    # Formateo a 2 decimais para caudais e 4 decimais para fraccións molares
    formatted_df = df_mesh_table.style.format({
        "Caudal Líquido L (kg/h)": "{:.2f}",
        "Caudal Vapor V (kg/h)": "{:.2f}",
        "Frac. Molar DME (x_DME)": "{:.4f}",
        "Frac. Molar MeOH (x_MeOH)": "{:.4f}",
        "Frac. Molar Auga (x_H2O)": "{:.4f}"
    })
    st.dataframe(formatted_df, use_container_width=True)

if btn_run:
    time.sleep(1.0)
    st.rerun()
