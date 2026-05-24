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
# INICIALIZACIÓN DO ESTADO
# ==========================================
if 'engine' not in st.session_state:
    st.session_state.engine = DigitalTwinEngine()
if 'analytics' not in st.session_state:
    st.session_state.analytics = ProcessAnalytics()
if 'is_running' not in st.session_state:
    st.session_state.is_running = False

# Variables de simulación P&ID
if 'pid_sim' not in st.session_state:
    st.session_state.pid_sim = {
        'time': [],
        'T_gasifier': [],
        'R_recycle': []
    }

engine = st.session_state.engine
analytics = st.session_state.analytics

def toggle_running():
    st.session_state.is_running = not st.session_state.is_running
    if st.session_state.is_running:
        st.session_state.caudal_dme = 438.47
        st.session_state.hotspot_t = 267.0
        st.session_state.presion_sistema = 110.0
    else:
        st.session_state.caudal_dme = 0.00
        st.session_state.hotspot_t = 150.0
        st.session_state.presion_sistema = 1.0

# ==========================================
# BARRA LATERAL (SIDEBAR)
# ==========================================
st.sidebar.title("🎛️ DCS Control")
st.sidebar.markdown("---")

if st.session_state.is_running:
    st.sidebar.button("⏸️ PAUSAR PLANTA", on_click=toggle_running, type="primary", use_container_width=True)
else:
    st.sidebar.button("▶️ ARRINCAR PLANTA", on_click=toggle_running, type="secondary", use_container_width=True)

if st.sidebar.button("🔄 RESET", use_container_width=True):
    st.session_state.engine = DigitalTwinEngine()
    st.session_state.is_running = False
    st.session_state.caudal_dme = 0.00
    st.session_state.hotspot_t = 150.0
    st.session_state.presion_sistema = 1.0
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
if abs(sum_comp - 100.0) > 0.1:
    st.sidebar.warning(f"Suma: {sum_comp:.2f}%. Debe ser 100%.")

if st.sidebar.button("💾 Aplicar Feedstock", use_container_width=True):
    engine.feedstock.set_custom_composition(c_val/100.0, h_val/100.0, o_val/100.0, n_val/100.0, s_val/100.0, ash_val/100.0)
    st.sidebar.success("Composición rexistrada.")

# Lóxica dinámica de simulación
if st.session_state.is_running:
    engine.simulation_step()
    
    # Alimentar o histórico do P&ID
    current_time = engine.t
    
    # Simulación de ruído/estabilización:
    # Setpoint T_gasifier = 750, comeza arredor e estabiliza
    t_gas_noise = random.uniform(-1.5, 1.5)
    t_gas = 750.0 + t_gas_noise
    
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

kpi_dme_delta = "Activo" if st.session_state.is_running else "Parado"

col1, col2, col3 = st.columns(3)
col1.metric("Caudal de DME (D3-1 Cabeza)", f"{st.session_state.caudal_dme:.2f} kg/h", delta=kpi_dme_delta)
col2.metric("Hotspot PFR Metanol (R2-1)", f"{st.session_state.hotspot_t:.1f} °C", delta=f"{st.session_state.hotspot_t - 150.0:.1f} °C vs base" if st.session_state.hotspot_t > 150 else "")
col3.metric("Descarga do Compresor (K2-1)", f"{st.session_state.presion_sistema:.1f} bar", delta=f"{st.session_state.presion_sistema - 1.0:.1f} bar vs Atm" if st.session_state.presion_sistema > 1 else "")

st.markdown("---")

# ==========================================
# PESTANAS DE VISUALIZACIÓN
# ==========================================
tab_scada, tab_pid, tab_analytics = st.tabs([
    "📊 SCADA Overview", 
    "🎛️ Control P&ID (Vivo)", 
    "📈 Motor Analítico e Balances"
])

# ----------------- TAB 1: SCADA -----------------
with tab_scada:
    col_left, col_right = st.columns(2)
    
    with col_left:
        st.markdown("### 🔥 Área Quente (Atm)")
        st.markdown(f"""
        **Sección de Gasificación**
        - Temperatura Reactor (RGibbs): `750.0 °C`
        - Presión de operación: `1.0 bar`
        - Caudal de Syngas Limpo (F1-1): `3057.81 kg/h`
        
        **Sección de Compresión (K2-1)**
        - Presión de Aspiración: `1.0 bar`
        - Presión de Descarga: `{st.session_state.presion_sistema:.1f} bar`
        - Etapas Activas: `5` (Ratio ~2.55)
        - Velocidade (RPM): `{'2500' if st.session_state.is_running else '0'}`
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
        - Pureza DME Cabeza D3-1: `99.90 %`
        """)

# ----------------- TAB 2: P&ID -----------------
with tab_pid:
    st.markdown("### Lazos de Control Dinámico (Anti-Windup)")
    
    if st.session_state.is_running and len(st.session_state.pid_sim['time']) > 0:
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
        st.dataframe(df_mass, use_container_width=True, hide_index=True)
        
        total_in = df_mass[df_mass["Tipo"] == "Entrada"]["Caudal (kg/h)"].sum()
        total_out = df_mass[df_mass["Tipo"] == "Saída"]["Caudal (kg/h)"].sum() + df_mass[df_mass["Tipo"] == "Inventario"]["Caudal (kg/h)"].sum()
        desviacion = abs(total_in - total_out)
        desv_pct = (desviacion / total_in) * 100.0 if total_in > 0 else 0.0
        st.success(f"✔️ **Balance pechado.** Desviación: {desviacion:.4f} kg/h ({desv_pct:.2f}%)")
        
    with col_table_right:
        st.markdown("#### Balance Elemental por Especies")
        st.dataframe(df_elem, use_container_width=True, hide_index=True)

if st.session_state.is_running:
    time.sleep(1.0)
    st.rerun()
