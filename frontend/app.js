// VARIABLES GLOBALES
let isRunning = false;
let simIntervalId = null;
const apiBase = ""; // Ruta relativa (funciona directamente servido polo backend)

// ELEMENTOS DOM
const btnPlay = document.getElementById("btn-play");
const btnPause = document.getElementById("btn-pause");
const btnStepMulti = document.getElementById("btn-step-multi");
const btnReset = document.getElementById("btn-reset");
const btnApply = document.getElementById("btn-apply-settings");

const txtTime = document.getElementById("txt-time");
const txtPressure = document.getElementById("txt-pressure");
const txtTempMax = document.getElementById("txt-temp-max");
const txtDmeFlow = document.getElementById("txt-dme-flow");
const txtValve = document.getElementById("txt-valve");
const txtFeedFlow = document.getElementById("txt-feed-flow");

const selectMode = document.getElementById("select-mode");
const sliderSetpoint = document.getElementById("slider-setpoint");
const valSetpoint = document.getElementById("val-setpoint");
const sliderCoolant = document.getElementById("slider-coolant");
const valCoolant = document.getElementById("val-coolant");
const sliderDisturbance = document.getElementById("slider-disturbance");
const valDisturbance = document.getElementById("val-disturbance");

const grpSetpoint = document.getElementById("grp-setpoint");
const grpCoolant = document.getElementById("grp-coolant");
const widgetTemp = document.getElementById("widget-temp");

// REFERENCIAS DE GRÁFICOS
let chartTimeseries = null;
let chartProfileTemp = null;
let chartProfileY = null;

// FUNCIÓNS AUXILIARES
function celsiusToKelvin(c) { return c + 273.15; }
function kelvinToCelsius(k) { return k - 273.15; }

// ACTUALIZACIÓN DE INTERFACE SEGUNDO MODO DE CONTROL
function updateUIForControlMode() {
    const mode = selectMode.value;
    if (mode === "MANUAL") {
        grpSetpoint.classList.add("disabled");
        sliderSetpoint.disabled = true;
        
        grpCoolant.classList.remove("disabled");
        sliderCoolant.disabled = false;
        sliderCoolant.classList.remove("text-muted");
        document.getElementById("val-coolant").classList.remove("text-muted");
    } else {
        // PID ou DMC
        grpSetpoint.classList.remove("disabled");
        sliderSetpoint.disabled = false;
        
        grpCoolant.classList.add("disabled");
        sliderCoolant.disabled = true;
        sliderCoolant.classList.add("text-muted");
        document.getElementById("val-coolant").classList.add("text-muted");
    }
}

// INICIALIZACIÓN DOS GRÁFICOS (CHART.JS)
function initCharts() {
    // Configuración global de fontes
    Chart.defaults.font.family = "'Outfit', sans-serif";
    Chart.defaults.color = '#a0a8c0';
    
    // 1. Gráfico de Evolución Temporal
    const ctxTs = document.getElementById("chart-timeseries").getContext("2d");
    chartTimeseries = new Chart(ctxTs, {
        type: 'line',
        data: {
            labels: [],
            datasets: [
                {
                    label: 'Tª Máx Reactor (°C)',
                    data: [],
                    borderColor: '#ff4d4d',
                    backgroundColor: 'rgba(255, 77, 77, 0.1)',
                    yAxisID: 'yTemp',
                    borderWidth: 2.5,
                    tension: 0.1,
                    pointRadius: 0
                },
                {
                    label: 'Tª Refrixerante (°C)',
                    data: [],
                    borderColor: '#00f2fe',
                    borderDash: [5, 5],
                    yAxisID: 'yTemp',
                    borderWidth: 2,
                    tension: 0.1,
                    pointRadius: 0
                },
                {
                    label: 'Setpoint Tª Máx (°C)',
                    data: [],
                    borderColor: '#ff9f1a',
                    borderDash: [3, 3],
                    yAxisID: 'yTemp',
                    borderWidth: 1.5,
                    tension: 0,
                    pointRadius: 0
                },
                {
                    label: 'Produción DME (mol/s)',
                    data: [],
                    borderColor: '#05c46b',
                    backgroundColor: 'rgba(5, 196, 107, 0.05)',
                    yAxisID: 'yFlow',
                    borderWidth: 2.5,
                    tension: 0.1,
                    pointRadius: 0,
                    fill: true
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            interaction: {
                mode: 'index',
                intersect: false,
            },
            plugins: {
                legend: {
                    position: 'top',
                    labels: { boxWidth: 15, padding: 15 }
                }
            },
            scales: {
                x: {
                    grid: { color: 'rgba(255, 255, 255, 0.03)' },
                    title: { display: true, text: 'Tempo (segundos)' }
                },
                yTemp: {
                    type: 'linear',
                    position: 'left',
                    grid: { color: 'rgba(255, 255, 255, 0.05)' },
                    title: { display: true, text: 'Temperatura (°C)' },
                    min: 220,
                    max: 290
                },
                yFlow: {
                    type: 'linear',
                    position: 'right',
                    grid: { drawOnChartArea: false },
                    title: { display: true, text: 'Caudal DME (mol/s)' },
                    min: 0,
                    max: 1.5
                }
            }
        }
    });

    // 2. Gráfico Perfil Axial de Temperatura
    const ctxProfTemp = document.getElementById("chart-profile-temp").getContext("2d");
    chartProfileTemp = new Chart(ctxProfTemp, {
        type: 'line',
        data: {
            labels: [],
            datasets: [{
                label: 'Temperatura (°C)',
                data: [],
                borderColor: '#ff9f1a',
                backgroundColor: 'rgba(255, 159, 26, 0.1)',
                borderWidth: 2.5,
                fill: true,
                tension: 0.2,
                pointBackgroundColor: '#ff9f1a'
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: { legend: { display: false } },
            scales: {
                x: {
                    grid: { color: 'rgba(255, 255, 255, 0.03)' },
                    title: { display: true, text: 'Posición axial (z, metros)' }
                },
                y: {
                    grid: { color: 'rgba(255, 255, 255, 0.05)' },
                    title: { display: true, text: 'Temperatura (°C)' },
                    min: 220,
                    max: 290
                }
            }
        }
    });

    // 3. Gráfico Perfil Axial de Composición
    const ctxProfY = document.getElementById("chart-profile-y").getContext("2d");
    chartProfileY = new Chart(ctxProfY, {
        type: 'line',
        data: {
            labels: [],
            datasets: [
                { label: 'H2', data: [], borderColor: '#4facfe', borderWidth: 2, tension: 0.1, pointRadius: 2 },
                { label: 'CO', data: [], borderColor: '#a0a8c0', borderWidth: 2, tension: 0.1, pointRadius: 2 },
                { label: 'CO2', data: [], borderColor: '#ff9f1a', borderWidth: 1.5, tension: 0.1, pointRadius: 0 },
                { label: 'MeOH', data: [], borderColor: '#e056fd', borderWidth: 2, tension: 0.1, pointRadius: 2 },
                { label: 'H2O', data: [], borderColor: '#22a6b3', borderWidth: 1.5, tension: 0.1, pointRadius: 0 },
                { label: 'DME', data: [], borderColor: '#05c46b', borderWidth: 3, tension: 0.1, pointRadius: 3 }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                legend: {
                    position: 'top',
                    labels: { boxWidth: 10, padding: 8 }
                }
            },
            scales: {
                x: {
                    grid: { color: 'rgba(255, 255, 255, 0.03)' },
                    title: { display: true, text: 'Posición axial (z, metros)' }
                },
                y: {
                    grid: { color: 'rgba(255, 255, 255, 0.05)' },
                    title: { display: true, text: 'Fracción molar (y_i)' },
                    min: 0,
                    max: 0.75
                }
            }
        }
    });
}

// ACTUALIZACIÓN DE GRÁFICOS E INDICADORES CON DATOS DA API
function updateDashboard(data) {
    // 1. Indicadores numéricos
    txtTime.innerText = `${data.time.toFixed(1)} s`;
    txtPressure.innerText = `${data.current_P_outlet.toFixed(2)} bar`;
    txtTempMax.innerText = `${data.current_T_max.toFixed(1)} °C`;
    txtDmeFlow.innerText = `${data.current_F_DME.toFixed(3)} mol/s`;
    txtValve.innerText = `${data.current_u_valve.toFixed(1)} %`;
    txtFeedFlow.innerText = `${data.current_F_feed.toFixed(2)} mol/s`;
    
    // Alarma por sobrequecemento / Runaway risco
    if (data.current_T_max > 275.0) {
        widgetTemp.classList.add("runaway-alert");
    } else {
        widgetTemp.classList.remove("runaway-alert");
    }
    
    // 2. Gráfico de series temporais
    chartTimeseries.data.labels = data.history.time.map(t => t.toFixed(0));
    chartTimeseries.data.datasets[0].data = data.history.T_max;
    chartTimeseries.data.datasets[1].data = data.history.T_coolant;
    
    // Setpoint replicado para a gráfica
    const spVal = data.T_setpoint;
    chartTimeseries.data.datasets[2].data = data.history.time.map((_, idx) => {
        // Se no momento histórico o modo era MANUAL, non debuxar setpoint
        return data.history.control_mode[idx] === "MANUAL" ? null : spVal;
    });
    
    chartTimeseries.data.datasets[3].data = data.history.F_DME_outlet;
    chartTimeseries.update('none'); // Update rápido sen animación para fluidez

    // 3. Perfís axiais
    const zLabels = data.profile.z.map(z => z.toFixed(2));
    
    chartProfileTemp.data.labels = zLabels;
    chartProfileTemp.data.datasets[0].data = data.profile.T;
    chartProfileTemp.update('none');

    chartProfileY.data.labels = zLabels;
    chartProfileY.data.datasets[0].data = data.profile.y_H2;
    chartProfileY.data.datasets[1].data = data.profile.y_CO;
    chartProfileY.data.datasets[2].data = data.profile.y_CO2;
    chartProfileY.data.datasets[3].data = data.profile.y_MeOH;
    chartProfileY.data.datasets[4].data = data.profile.y_H2O;
    chartProfileY.data.datasets[5].data = data.profile.y_DME;
    chartProfileY.update('none');
}

// PETICIÓNS API
async function fetchStatus() {
    try {
        const response = await fetch(`${apiBase}/api/status`);
        if (!response.ok) throw new Error("Erro ao obter estado");
        const data = await response.json();
        
        // Cando cargamos o estado por primeira vez, poñer os valores nos controis
        if (!isRunning) {
            selectMode.value = data.control_mode;
            sliderSetpoint.value = data.T_setpoint;
            valSetpoint.innerText = `${data.T_setpoint.toFixed(1)} °C`;
            
            sliderCoolant.value = data.current_T_coolant;
            valCoolant.innerText = `${data.current_T_coolant.toFixed(1)} °C`;
            
            sliderDisturbance.value = data.disturbance_flow_pct;
            valDisturbance.innerText = `${data.disturbance_flow_pct > 0 ? '+' : ''}${data.disturbance_flow_pct} %`;
            
            updateUIForControlMode();
        }
        
        updateDashboard(data);
    } catch (err) {
        console.error(err);
    }
}

async function sendSettings() {
    const payload = {
        control_mode: selectMode.value,
        T_setpoint: parseFloat(sliderSetpoint.value),
        T_coolant_manual: parseFloat(sliderCoolant.value),
        disturbance_flow_pct: parseFloat(sliderDisturbance.value)
    };
    
    try {
        const response = await fetch(`${apiBase}/api/settings`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });
        if (!response.ok) throw new Error("Erro ao gardar axustes");
        
        // Solicitar estado actualizado
        await fetchStatus();
    } catch (err) {
        console.error("Erro ao aplicar axustes: ", err);
    }
}

async function simStep() {
    try {
        const response = await fetch(`${apiBase}/api/step`, { method: 'POST' });
        if (!response.ok) throw new Error("Erro no paso de simulación");
        await fetchStatus();
    } catch (err) {
        console.error(err);
        pauseSimulation();
    }
}

async function simStepMulti() {
    try {
        // Bloquear temporalmente o botón
        btnStepMulti.disabled = true;
        const response = await fetch(`${apiBase}/api/step-multi?steps=5`, { method: 'POST' });
        if (!response.ok) throw new Error("Erro no paso múltiple");
        await fetchStatus();
    } catch (err) {
        console.error(err);
    } finally {
        btnStepMulti.disabled = false;
    }
}

async function resetSimulation() {
    if (confirm("Está seguro de que quere reiniciar a simulación ao estado térmico inicial?")) {
        pauseSimulation();
        try {
            const response = await fetch(`${apiBase}/api/reset`, { method: 'POST' });
            if (!response.ok) throw new Error("Erro ao reiniciar");
            await fetchStatus();
        } catch (err) {
            console.error(err);
        }
    }
}

// XESTIÓN DE LAZO DE SIMULACIÓN
function playSimulation() {
    if (isRunning) return;
    isRunning = true;
    btnPlay.disabled = true;
    btnPause.disabled = false;
    
    // Lazo de paso temporal constante (cada 1s simúlase Ts=2s do xemelo)
    simIntervalId = setInterval(simStep, 1000);
}

function pauseSimulation() {
    if (!isRunning) return;
    isRunning = false;
    btnPlay.disabled = false;
    btnPause.disabled = true;
    
    clearInterval(simIntervalId);
    simIntervalId = null;
}

// ASOCIACIÓN DE EVENTOS
selectMode.addEventListener("change", () => {
    updateUIForControlMode();
    sendSettings();
});

sliderSetpoint.addEventListener("input", (e) => {
    valSetpoint.innerText = `${parseFloat(e.target.value).toFixed(1)} °C`;
});
sliderSetpoint.addEventListener("change", sendSettings);

sliderCoolant.addEventListener("input", (e) => {
    valCoolant.innerText = `${parseFloat(e.target.value).toFixed(1)} °C`;
});
sliderCoolant.addEventListener("change", sendSettings);

sliderDisturbance.addEventListener("input", (e) => {
    const val = parseInt(e.target.value);
    valDisturbance.innerText = `${val > 0 ? '+' : ''}${val} %`;
});
sliderDisturbance.addEventListener("change", sendSettings);

btnPlay.addEventListener("click", playSimulation);
btnPause.addEventListener("click", pauseSimulation);
btnStepMulti.addEventListener("click", simStepMulti);
btnReset.addEventListener("click", resetSimulation);
btnApply.addEventListener("click", sendSettings);

// EVENTO DE CARGA INICIAL
window.addEventListener("DOMContentLoaded", () => {
    initCharts();
    fetchStatus();
    updateUIForControlMode();
});
