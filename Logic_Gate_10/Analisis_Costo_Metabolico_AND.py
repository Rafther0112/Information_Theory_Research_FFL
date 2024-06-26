#%%
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

Red_color = "#d62728"
Orange_color = "#ff7f0e"
Green_color = "#2ca02c"
Blue_color =  "#1f77b4"

# %% Cargamos los datos de las simulaciones para usarlos con OR
simulacion_FFL_C1 = np.load('Resultados_Simulacion/Simulacion_FFL_C1_AND_final.npy', allow_pickle=True).item()
simulacion_FFL_C2 = np.load('Resultados_Simulacion/Simulacion_FFL_C2_AND_final.npy', allow_pickle=True).item()
simulacion_FFL_C3 = np.load('Resultados_Simulacion/Simulacion_FFL_C3_AND_final.npy', allow_pickle=True).item()
simulacion_FFL_C4 = np.load('Resultados_Simulacion/Simulacion_FFL_C4_AND_final.npy', allow_pickle=True).item()

simulacion_FFL_I1 = np.load('Resultados_Simulacion/Simulacion_FFL_I1_AND_final.npy', allow_pickle=True).item()
simulacion_FFL_I2 = np.load('Resultados_Simulacion/Simulacion_FFL_I2_AND_final.npy', allow_pickle=True).item()
simulacion_FFL_I3 = np.load('Resultados_Simulacion/Simulacion_FFL_I3_AND_final.npy', allow_pickle=True).item()
simulacion_FFL_I4 = np.load('Resultados_Simulacion/Simulacion_FFL_I4_AND_final.npy', allow_pickle=True).item()
#%%
# _________________ANALISIS DE COSTO DE PROTEINAS DE FUNCIONAMIENTO____________________________
#%% Hacemos la lectura de cantidad de proteina en estado estacionario coherentes
sampling = 150
Costo_Metabolico_C1_AND = np.empty(shape=(4,10), dtype='float')
Costo_Metabolico_C2_AND = np.empty(shape=(4,10), dtype='float')
Costo_Metabolico_C3_AND = np.empty(shape=(4,10), dtype='float')
Costo_Metabolico_C4_AND = np.empty(shape=(4,10), dtype='float')

for posicion_Hill, Valor_Hill in enumerate(range(1,5)):
    for posicion_Valor_K, valor_K in enumerate(range(0,10)):
        Costo_Metabolico_C1_AND[posicion_Hill][posicion_Valor_K] = (np.mean(np.mean(simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]) + np.mean(np.mean(simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))/195000
        Costo_Metabolico_C2_AND[posicion_Hill][posicion_Valor_K] = (np.mean(np.mean(simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]) + np.mean(np.mean(simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))/195000
        Costo_Metabolico_C3_AND[posicion_Hill][posicion_Valor_K] = (np.mean(np.mean(simulacion_FFL_C3[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]) + np.mean(np.mean(simulacion_FFL_C3[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))/195000
        Costo_Metabolico_C4_AND[posicion_Hill][posicion_Valor_K] = (np.mean(np.mean(simulacion_FFL_C4[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]) + np.mean(np.mean(simulacion_FFL_C4[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))/195000


np.save('Resultados_Costos_Metabolicos/Costo_Metabolico_C1_AND.npy', Costo_Metabolico_C1_AND)
np.save('Resultados_Costos_Metabolicos/Costo_Metabolico_C2_AND.npy', Costo_Metabolico_C2_AND)
np.save('Resultados_Costos_Metabolicos/Costo_Metabolico_C3_AND.npy', Costo_Metabolico_C3_AND)
np.save('Resultados_Costos_Metabolicos/Costo_Metabolico_C4_AND.npy', Costo_Metabolico_C4_AND)

#%%
#%% Hacemos la lectura de cantidad de proteina en estado estacionario coherentes
sampling = 150
Costo_Metabolico_I1_AND = np.empty(shape=(4,10), dtype='float')
Costo_Metabolico_I2_AND = np.empty(shape=(4,10), dtype='float')
Costo_Metabolico_I3_AND = np.empty(shape=(4,10), dtype='float')
Costo_Metabolico_I4_AND = np.empty(shape=(4,10), dtype='float')

for posicion_Hill, Valor_Hill in enumerate(range(1,5)):
    for posicion_Valor_K, valor_K in enumerate(range(0,10)):
        Costo_Metabolico_I1_AND[posicion_Hill][posicion_Valor_K] = (np.mean(np.mean(simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]) + np.mean(np.mean(simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))/195000
        Costo_Metabolico_I2_AND[posicion_Hill][posicion_Valor_K] = (np.mean(np.mean(simulacion_FFL_I2[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]) + np.mean(np.mean(simulacion_FFL_I2[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))/195000
        Costo_Metabolico_I3_AND[posicion_Hill][posicion_Valor_K] = (np.mean(np.mean(simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]) + np.mean(np.mean(simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))/195000
        Costo_Metabolico_I4_AND[posicion_Hill][posicion_Valor_K] = (np.mean(np.mean(simulacion_FFL_I4[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]) + np.mean(np.mean(simulacion_FFL_I4[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))/195000

np.save('Resultados_Costos_Metabolicos/Costo_Metabolico_I1_AND.npy', Costo_Metabolico_I1_AND)
np.save('Resultados_Costos_Metabolicos/Costo_Metabolico_I2_AND.npy', Costo_Metabolico_I2_AND)
np.save('Resultados_Costos_Metabolicos/Costo_Metabolico_I3_AND.npy', Costo_Metabolico_I3_AND)
np.save('Resultados_Costos_Metabolicos/Costo_Metabolico_I4_AND.npy', Costo_Metabolico_I4_AND)
# %%
