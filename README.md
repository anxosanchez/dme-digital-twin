# Xemelo Dixital: Planta de Síntese de DME a partir de Biomasa

Este repositorio contén o motor de simulación e o panel de control (Dashboard) dun Xemelo Dixital dinámico desenvolvido en Python, para modelar e controlar unha planta de produción de Éter Dimetílico (DME).

O modelo matemático foi calibrado integramente coas simulacións e as propiedades físicas proporcionadas por Aspen Plus, incorporando:
- **Gasificación de Biomasa:** Modelo de equilibrio termodinámico (Minimización da Enerxía Libre de Gibbs).
- **Reactor PFR de Metanol:** Cinética heteroxénea estrita LHHW para as reaccións de síntese de metanol e WGS inversa.
- **Torres de Destilación e Separación:** Modelos de separación termodinámica con ecuacións de estado (Peng-Robinson EoS).
- **Control Automático Avanzado:** Controladores PID e DMC (Predictivo) integrados nativamente para corrixir os reciclos, acumulación de inertes (Nitróxeno) e xestión térmica da planta.

## Estrutura do Proxecto
- `app.py`: A interface interactiva SCADA baseada en Streamlit. Actúa como o cerebro central do panel.
- `engine.py`: O "Engine" do xemelo dixital, que executa o bucle temporal dinámico conectando tódolos equipos.
- `analytics.py`: Módulo adicado á xeración dos balances de materia, gráficos térmicos e perfís axiais do reactor.
- E os módulos base para cada operación unitaria: `gasification.py`, `reactor.py`, `methanol_reactor.py`, `dme_reactor.py`, `gas_compression.py`, `separation.py` e `control.py`.

## Como executar o Proxecto

Asegúrate de ter Python instalado no teu sistema. Para executar o xemelo dixital en local:

1. Instala as dependencias do proxecto:
   ```bash
   pip install -r requirements.txt
   ```

2. Arranca o servidor de Streamlit:
   ```bash
   streamlit run app.py
   ```

3. O teu navegador debería abrir automaticamente a ruta `http://localhost:8501`. Dende aí, introduce os valores da biomasa e fai clic en **Arrincar Planta**.
