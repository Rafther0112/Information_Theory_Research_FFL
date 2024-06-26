#%%
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

Red_color = "#d62728"
Orange_color = "#ff7f0e"
Green_color = "#2ca02c"
Blue_color =  "#1f77b4"

def conditional_covarianza_Z_Y(covariance_matrix):
    return covariance_matrix[2][2] - ((covariance_matrix[2][1])**2)/(covariance_matrix[1][1])

def conditional_covarianza_X_ZdadoY(covariance_matrix):
    return covariance_matrix[0][2] - ((covariance_matrix[0][1])*(covariance_matrix[2][1]))/(covariance_matrix[1][1])

# %% Cargamos los datos de las simulaciones para usarlos con OR
simulacion_FFL_C1 = np.load('Resultados_Simulacion/Simulacion_FFL_C1_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_C2 = np.load('Resultados_Simulacion/Simulacion_FFL_C2_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_C3 = np.load('Resultados_Simulacion/Simulacion_FFL_C3_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_C4 = np.load('Resultados_Simulacion/Simulacion_FFL_C4_OR_final.npy', allow_pickle=True).item()

simulacion_FFL_I1 = np.load('Resultados_Simulacion/Simulacion_FFL_I1_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_I2 = np.load('Resultados_Simulacion/Simulacion_FFL_I2_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_I3 = np.load('Resultados_Simulacion/Simulacion_FFL_I3_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_I4 = np.load('Resultados_Simulacion/Simulacion_FFL_I4_OR_final.npy', allow_pickle=True).item()
#%%
# _________________ANALISIS DE COSTO DE PROTEINAS DE FUNCIONAMIENTO____________________________
#%% Hacemos la lectura de cantidad de proteina en estado estacionario coherentes
sampling = 150
proteinas_estado_estacionario_Y_C = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]
proteinas_estado_estacionario_Z_C = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]

for posicion_Hill, Valor_Hill in enumerate(range(1,5)):
    for valor_K in range(0,10):
        proteinas_estado_estacionario_Y_C[posicion_Hill][0].append(np.mean(np.mean(simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Y_C[posicion_Hill][1].append(np.mean(np.mean(simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Y_C[posicion_Hill][2].append(np.mean(np.mean(simulacion_FFL_C3[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Y_C[posicion_Hill][3].append(np.mean(np.mean(simulacion_FFL_C4[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]))

        proteinas_estado_estacionario_Z_C[posicion_Hill][0].append(np.mean(np.mean(simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Z_C[posicion_Hill][1].append(np.mean(np.mean(simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Z_C[posicion_Hill][2].append(np.mean(np.mean(simulacion_FFL_C3[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Z_C[posicion_Hill][3].append(np.mean(np.mean(simulacion_FFL_C4[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))
#%% Grafia para los coherentes
regulacion_K = 0
plt.figure(figsize=(8,6))
plt.title(r"Costo metabólico FFL Coherentes general" + "\n" + r"compuerta lógica OR", fontsize = 16)
plt.xlabel(r"Valor tasa $K_x$", fontsize = 16)
plt.ylabel(r" $\langle Y \rangle + \langle Z \rangle$", fontsize = 16)

plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), label = "C1", color = Red_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), label = "C2", color = Orange_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), label = "C3", color = Green_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), label = "C4", color = Blue_color)

regulacion_K = 1
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), color = Red_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), color = Orange_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), color = Green_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), color = Blue_color)

regulacion_K = 2
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), color = Red_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), color = Orange_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), color = Green_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), color = Blue_color)

regulacion_K = 3
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), color = Red_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), color = Orange_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), color = Green_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), color = Blue_color)

plt.axhline(y = 195000, linestyle = "--", color = "gray", label = "Hill = 0")
plt.legend()
#plt.savefig("Imagenes_Resultados/Costo_Metabolico_Coherente_General_OR.jpg", dpi=1000)
#%% Graficas multiples
# Crea una figura y un arreglo de ejes de 2x2
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
fig.suptitle(r"Costo metabólico FFL Coherentes con diferente coeficiente de Hill" + "\n" + r"compuerta lógica OR", fontsize=16)

regulacion_K = 0
# Puedes acceder a cada eje individualmente utilizando la notación de índices
axes[0, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), label = "C1", color = Red_color, alpha = 0.4)
axes[0, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), label = "C2", color = Orange_color, alpha = 0.4)
axes[0, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), label = "C3", color = Green_color, alpha = 0.4)
axes[0, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), label = "C4", color = Blue_color, alpha = 0.4)
axes[0, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), color = Red_color, linestyle = "dashed")
axes[0, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), color = Orange_color, linestyle = "dashed")
axes[0, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), color = Green_color, linestyle = "dashed")
axes[0, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), color = Blue_color, linestyle = "dashed")

axes[0, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 0].set_ylabel(r" $\langle Y \rangle + \langle Z \rangle$", fontsize = 14)
axes[0, 0].set_ylim([189000, 200000]) 
axes[0, 0].axhline(y = 195000, linestyle = "--", color = "gray", label = "Hill = 0")
axes[0, 0].legend() 

regulacion_K = 1
axes[0, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), label = "C1", color = Red_color, alpha = 0.4)
axes[0, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), label = "C2", color = Orange_color, alpha = 0.4)
axes[0, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), label = "C3", color = Green_color, alpha = 0.4)
axes[0, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), label = "C4", color = Blue_color, alpha = 0.4)
axes[0, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), color = Red_color, linestyle = "dashed")
axes[0, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), color = Orange_color, linestyle = "dashed")
axes[0, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), color = Green_color, linestyle = "dashed")
axes[0, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), color = Blue_color, linestyle = "dashed")

axes[0, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 1].set_ylabel(r" $\langle Y \rangle + \langle Z \rangle$", fontsize = 14)
axes[0, 1].set_ylim([189000, 200000])
axes[0, 1].axhline(y = 195000, linestyle = "--", color = "gray", label = "Hill = 0")  
axes[0, 1].legend()

regulacion_K = 2
axes[1, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), label = "C1", color = Red_color, alpha = 0.4)
axes[1, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), label = "C2", color = Orange_color, alpha = 0.4)
axes[1, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), label = "C3", color = Green_color, alpha = 0.4)
axes[1, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), label = "C4", color = Blue_color, alpha = 0.4)
axes[1, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), color = Red_color, linestyle = "dashed")
axes[1, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), color = Orange_color, linestyle = "dashed")
axes[1, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), color = Green_color, linestyle = "dashed")
axes[1, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), color = Blue_color, linestyle = "dashed")

axes[1, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 0].set_ylabel(r" $\langle Y \rangle + \langle Z \rangle$", fontsize = 14)
axes[1, 0].set_ylim([189000, 204500])  
axes[1, 0].axhline(y = 195000, linestyle = "--", color = "gray", label = "Hill = 0" )
axes[1, 0].legend()

regulacion_K = 3
axes[1, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), label = "C1", color = Red_color, alpha = 0.4)
axes[1, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), label = "C2", color = Orange_color, alpha = 0.4)
axes[1, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), label = "C3", color = Green_color, alpha = 0.4)
axes[1, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), label = "C4", color = Blue_color, alpha = 0.4)
axes[1, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][0]), color = Red_color, linestyle = "dashed")
axes[1, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][1]), color = Orange_color, linestyle = "dashed")
axes[1, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][2]), color = Green_color, linestyle = "dashed")
axes[1, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_C[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_C[regulacion_K][3]), color = Blue_color, linestyle = "dashed")

axes[1, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 1].set_ylabel(r" $\langle Y \rangle + \langle Z \rangle$", fontsize = 14)
axes[1, 1].set_ylim([189000, 204500])  
axes[1, 1].axhline(y = 195000, linestyle = "--", color = "gray", label = "Hill = 0")
axes[1, 1].legend()

plt.tight_layout()
#plt.savefig("Imagenes_Resultados/Costo_Metabolico_Coherente_Multiple_OR.jpg", dpi=1000)
#%%

# %% Lectura de cantidad de proteina en estado estacionario Incoherentes
proteinas_estado_estacionario_Y_I = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]
proteinas_estado_estacionario_Z_I = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]
sampling = 150
for posicion_Hill, Valor_Hill in enumerate(range(1,5)):
    for valor_K in range(0,10):
        proteinas_estado_estacionario_Y_I[posicion_Hill][0].append(np.mean(np.mean(simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Y_I[posicion_Hill][1].append(np.mean(np.mean(simulacion_FFL_I2[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Y_I[posicion_Hill][2].append(np.mean(np.mean(simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Y_I[posicion_Hill][3].append(np.mean(np.mean(simulacion_FFL_I4[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K], axis=0)[sampling:-1]))

        proteinas_estado_estacionario_Z_I[posicion_Hill][0].append(np.mean(np.mean(simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Z_I[posicion_Hill][1].append(np.mean(np.mean(simulacion_FFL_I2[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Z_I[posicion_Hill][2].append(np.mean(np.mean(simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))
        proteinas_estado_estacionario_Z_I[posicion_Hill][3].append(np.mean(np.mean(simulacion_FFL_I4[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K], axis=0)[sampling:-1]))
# %% Grafica para los incoherentes

regulacion_K = 0
plt.figure(figsize=(8,6))
plt.title(r"Costo metabólico FFL Incoherentes general" + "\n" + r"compuerta lógica OR", fontsize = 16)
plt.xlabel(r"Valor tasa $K_x$", fontsize = 16)
plt.ylabel(r" $\langle Y \rangle + \langle Z \rangle$", fontsize = 16)

plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), label = "I1", color = Red_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), label = "I2", color = Orange_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), label = "I3", color = Green_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), label = "I4", color = Blue_color)

regulacion_K = 1
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), color = Red_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), color = Orange_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), color = Green_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), color = Blue_color)

regulacion_K = 2
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), color = Red_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), color = Orange_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), color = Green_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), color = Blue_color)

regulacion_K = 3
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), color = Red_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), color = Orange_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), color = Green_color)
plt.scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), color = Blue_color)

plt.axhline(y = 195000, linestyle = "--", color = "gray", label = "Hill = 0")
plt.legend()
plt.savefig("Imagenes_Resultados/Costo_Metabolico_Incoherente_General_OR.jpg", dpi=1000)
#%% Grafica multiple para los incoherentes
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
fig.suptitle(r"Costo metabólico FFL Incoherentes con diferente coeficiente de Hill" + "\n" + r"compuerta lógica OR", fontsize=16)

regulacion_K = 0
# Puedes acceder a cada eje individualmente utilizando la notación de índices
axes[0, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), label = "I1", color = Red_color, alpha = 0.4)
axes[0, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), label = "I2", color = Orange_color, alpha = 0.4)
axes[0, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), label = "I3", color = Green_color, alpha = 0.4)
axes[0, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), label = "I4", color = Blue_color, alpha = 0.4)
axes[0, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), linestyle = "dashed", color = Red_color)
axes[0, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), linestyle = "dashed", color = Orange_color)
axes[0, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), linestyle = "dashed", color = Green_color)
axes[0, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), linestyle = "dashed", color = Blue_color)

axes[0, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 0].set_ylabel(r" $\langle Y \rangle + \langle Z \rangle$", fontsize = 14)
axes[0, 0].set_ylim([192000, 200000]) 
axes[0, 0].axhline(y = 195000, linestyle = "--", color = "gray", label = "Hill = 0")
axes[0, 0].legend() 

regulacion_K = 1
axes[0, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), label = "I1", color = Red_color, alpha = 0.4)
axes[0, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), label = "I2", color = Orange_color, alpha = 0.4)
axes[0, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), label = "I3", color = Green_color, alpha = 0.4)
axes[0, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), label = "I4", color = Blue_color, alpha = 0.4)
axes[0, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), linestyle = "dashed", color = Red_color)
axes[0, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), linestyle = "dashed", color = Orange_color)
axes[0, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), linestyle = "dashed", color = Green_color)
axes[0, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), linestyle = "dashed", color = Blue_color)

axes[0, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 1].set_ylabel(r" $\langle Y \rangle + \langle Z \rangle$", fontsize = 14)
axes[0, 1].set_ylim([192000, 200000])
axes[0, 1].axhline(y = 195000, linestyle = "--", color = "gray", label = "Hill = 0")  
axes[0, 1].legend()

regulacion_K = 2
axes[1, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), label = "I1", color = Red_color, alpha = 0.4)
axes[1, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), label = "I2", color = Orange_color, alpha = 0.4)
axes[1, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), label = "I3", color = Green_color, alpha = 0.4)
axes[1, 0].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), label = "I4", color = Blue_color, alpha = 0.4)
axes[1, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), linestyle = "dashed", color = Red_color)
axes[1, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), linestyle = "dashed", color = Orange_color)
axes[1, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), linestyle = "dashed", color = Green_color)
axes[1, 0].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), linestyle = "dashed", color = Blue_color)

axes[1, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 0].set_ylabel(r" $\langle Y \rangle + \langle Z \rangle$", fontsize = 14)
axes[1, 0].set_ylim([189000, 200500])  
axes[1, 0].axhline(y = 195000, linestyle = "--", color = "gray", label = "Hill = 0" )
axes[1, 0].legend()

regulacion_K = 3
axes[1, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), label = "I1", color = Red_color, alpha = 0.4)
axes[1, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), label = "I2", color = Orange_color, alpha = 0.4)
axes[1, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), label = "I3", color = Green_color, alpha = 0.4)
axes[1, 1].scatter(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), label = "I4", color = Blue_color, alpha = 0.4)
axes[1, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][0]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][0]), linestyle = "dashed", color = Red_color)
axes[1, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][1]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][1]), linestyle = "dashed", color = Orange_color)
axes[1, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][2]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][2]), linestyle = "dashed", color = Green_color)
axes[1, 1].plot(range(1,11), np.array(proteinas_estado_estacionario_Y_I[regulacion_K][3]) + np.array(proteinas_estado_estacionario_Z_I[regulacion_K][3]), linestyle = "dashed", color = Blue_color)

axes[1, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 1].set_ylabel(r" $\langle Y \rangle + \langle Z \rangle$", fontsize = 14)
axes[1, 1].set_ylim([189000, 200500])  
axes[1, 1].axhline(y = 195000, linestyle = "--", color = "gray", label = "Hill = 0")
axes[1, 1].legend()
plt.tight_layout()
plt.savefig("Imagenes_Resultados/Costo_Metabolico_Incoherente_Multiple_OR.jpg", dpi=1000)
# %%
# _________________ANALISIS DE INFORMACION EN ESTADO ESTACIONARIO____________________________
#%% Analisis para los coherentes
sampling = 150
informacion_estado_estacionario_C_X_Z = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]
informacion_estado_estacionario_C_Z_X_Y = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]

for posicion_Hill, Valor_Hill in enumerate(range(1,5)):
    for valor_K in range(0,10):

        data_C1 = {'X': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C1 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C1)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][0].append((1/2)*np.log2((Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2])/(Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2] - (Cov_matrix_C1[0][2])**2)))
        informacion_estado_estacionario_C_Z_X_Y[posicion_Hill][0].append((1/2)*np.log2((Cov_matrix_C1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C1))/((Cov_matrix_C1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C1))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C1))**2)))

        data_C2 = {'X': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C2 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C2)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][1].append((1/2)*np.log2((Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2])/(Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2] - (Cov_matrix_C2[0][2])**2)))
        informacion_estado_estacionario_C_Z_X_Y[posicion_Hill][1].append((1/2)*np.log2((Cov_matrix_C2[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C2))/((Cov_matrix_C2[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C2))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C2))**2)))

        data_C3 = {'X': simulacion_FFL_C3[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_C3[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_C3[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C3 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C3)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][2].append((1/2)*np.log2((Cov_matrix_C3[0][0]* Cov_matrix_C3[2][2])/(Cov_matrix_C3[0][0]* Cov_matrix_C3[2][2] - (Cov_matrix_C3[0][2])**2)))
        informacion_estado_estacionario_C_Z_X_Y[posicion_Hill][2].append((1/2)*np.log2((Cov_matrix_C3[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C3))/((Cov_matrix_C3[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C3))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C3))**2)))

        data_C4 = {'X': simulacion_FFL_C4[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_C4[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_C4[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C4 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C4)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][3].append((1/2)*np.log2((Cov_matrix_C4[0][0]* Cov_matrix_C4[2][2])/(Cov_matrix_C4[0][0]* Cov_matrix_C4[2][2] - (Cov_matrix_C4[0][2])**2)))
        informacion_estado_estacionario_C_Z_X_Y[posicion_Hill][3].append((1/2)*np.log2((Cov_matrix_C4[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C4))/((Cov_matrix_C4[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C4))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C4))**2)))


# %% Graficas multiples de informacion coherentes I(X;Z)
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
fig.suptitle(r"Información mutua en estado estacionario $I(\langle X \rangle;\langle Z \rangle)$" + "\n" +  r"  FFL Coherentes con diferente coeficiente de Hill |  Compuerta lógica OR", fontsize = 16)

regulacion_K = 0
# Puedes acceder a cada eje individualmente utilizando la notación de índices
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "C1", color = Red_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "C2", color = Orange_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "C3", color = Green_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "C4", color = Blue_color)
axes[0, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[0, 0].set_ylim([-0.005, 0.09]) 
axes[0, 0].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[0, 0].legend() 

regulacion_K = 1
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "C1", color = Red_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "C2", color = Orange_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "C3", color = Green_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "C4", color = Blue_color)
axes[0, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[0, 1].set_ylim([-0.005, 0.09]) 
axes[0, 1].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[0, 1].legend()

regulacion_K = 2
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "C1", color = Red_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "C2", color = Orange_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "C3", color = Green_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "C4", color = Blue_color)
axes[1, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[1, 0].set_ylim([-0.005, 0.14]) 
axes[1, 0].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[1, 0].legend()

regulacion_K = 3
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "C1", color = Red_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "C2", color = Orange_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "C3", color = Green_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "C4", color = Blue_color)
axes[1, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[1, 1].set_ylim([-0.005, 0.14]) 
axes[1, 1].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[1, 1].legend()

plt.tight_layout()
#plt.savefig("Imagenes_Resultados/Informacion_Estacionaria_Coherente_Multiple_I(X_Z)_OR.jpg", dpi=1000)
# %% Graficas multiples de informacion coherentes I(X;Z|Y)

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
fig.suptitle(r"Información mutua en estado estacionario $I(\langle X \rangle;\langle Z \rangle | \langle Y \rangle)$" + "\n" + r"  FFL Coherentes con diferente coeficiente de Hill |  Compuerta lógica OR", fontsize = 16)

regulacion_K = 0
# Puedes acceder a cada eje individualmente utilizando la notación de índices
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][0]), label = "C1", color = Red_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][1]), label = "C2", color = Orange_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][2]), label = "C3", color = Green_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][3]), label = "C4", color = Blue_color)
axes[0, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle | \langle Y \rangle)$", fontsize = 14)
axes[0, 0].set_ylim([-0.001, 0.0175]) 
axes[0, 0].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[0, 0].legend() 

regulacion_K = 1
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][0]), label = "C1", color = Red_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][1]), label = "C2", color = Orange_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][2]), label = "C3", color = Green_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][3]), label = "C4", color = Blue_color)
axes[0, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle | \langle Y \rangle)$", fontsize = 14)
axes[0, 1].set_ylim([-0.001, 0.0175]) 
axes[0, 1].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[0, 1].legend()

regulacion_K = 2
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][0]), label = "C1", color = Red_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][1]), label = "C2", color = Orange_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][2]), label = "C3", color = Green_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][3]), label = "C4", color = Blue_color)
axes[1, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle | \langle Y \rangle)$", fontsize = 14)
axes[1, 0].set_ylim([-0.001, 0.0275]) 
axes[1, 0].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[1, 0].legend()

regulacion_K = 3
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][0]), label = "C1", color = Red_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][1]), label = "C2", color = Orange_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][2]), label = "C3", color = Green_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][3]), label = "C4", color = Blue_color)
axes[1, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle| \langle Y \rangle)$", fontsize = 14)
axes[1, 1].set_ylim([-0.001, 0.0275]) 
axes[1, 1].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[1, 1].legend()

plt.tight_layout()
plt.savefig("Imagenes_Resultados/Informacion_Estacionaria_Coherente_Multiple_I(X_Z|Y)_OR.jpg", dpi=1000)

#%% Analisis para los incoherentes
sampling = 100
informacion_estado_estacionario_C_X_Z = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]
informacion_estado_estacionario_C_Z_X_Y = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]

for posicion_Hill, Valor_Hill in enumerate(range(1,5)):
    for valor_K in range(0,10):

        data_C1 = {'X': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C1 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C1)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][0].append((1/2)*np.log2((Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2])/(Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2] - (Cov_matrix_C1[0][2])**2)))
        informacion_estado_estacionario_C_Z_X_Y[posicion_Hill][0].append((1/2)*np.log2((Cov_matrix_C1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C1))/((Cov_matrix_C1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C1))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C1))**2)))

        data_C2 = {'X': simulacion_FFL_I2[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_I2[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_I2[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C2 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C2)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][1].append((1/2)*np.log2((Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2])/(Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2] - (Cov_matrix_C2[0][2])**2)))
        informacion_estado_estacionario_C_Z_X_Y[posicion_Hill][1].append((1/2)*np.log2((Cov_matrix_C2[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C2))/((Cov_matrix_C2[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C2))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C2))**2)))

        data_C3 = {'X': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C3 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C3)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][2].append((1/2)*np.log2((Cov_matrix_C3[0][0]* Cov_matrix_C3[2][2])/(Cov_matrix_C3[0][0]* Cov_matrix_C3[2][2] - (Cov_matrix_C3[0][2])**2)))
        informacion_estado_estacionario_C_Z_X_Y[posicion_Hill][2].append((1/2)*np.log2((Cov_matrix_C3[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C3))/((Cov_matrix_C3[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C3))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C3))**2)))

        data_C4 = {'X': simulacion_FFL_I4[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_I4[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_I4[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C4 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C4)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][3].append((1/2)*np.log2((Cov_matrix_C4[0][0]* Cov_matrix_C4[2][2])/(Cov_matrix_C4[0][0]* Cov_matrix_C4[2][2] - (Cov_matrix_C4[0][2])**2)))
        informacion_estado_estacionario_C_Z_X_Y[posicion_Hill][3].append((1/2)*np.log2((Cov_matrix_C4[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C4))/((Cov_matrix_C4[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C4))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C4))**2)))
# %% Graficas multiples de informacion Incoherentes I(X;Z)

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
fig.suptitle(r"Información mutua en estado estacionario $I(\langle X \rangle;\langle Z \rangle)$" + "\n" +  r"  FFL Incoherentes con diferente coeficiente de Hill |  Compuerta lógica OR", fontsize = 16)

regulacion_K = 0
# Puedes acceder a cada eje individualmente utilizando la notación de índices
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "I1", color = Red_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "I2", color = Orange_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "I3", color = Green_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "I4", color = Blue_color)
axes[0, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[0, 0].set_ylim([-0.005, 0.08]) 
axes[0, 0].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[0, 0].legend() 

regulacion_K = 1
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "I1", color = Red_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "I2", color = Orange_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "I3", color = Green_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "I4", color = Blue_color)
axes[0, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[0, 1].set_ylim([-0.005, 0.08]) 
axes[0, 1].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[0, 1].legend()

regulacion_K = 2
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "I1", color = Red_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "I2", color = Orange_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "I3", color = Green_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "I4", color = Blue_color)
axes[1, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[1, 0].set_ylim([-0.005, 0.2]) 
axes[1, 0].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[1, 0].legend()

regulacion_K = 3
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "I1", color = Red_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "I2", color = Orange_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "I3", color = Green_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "I4", color = Blue_color)
axes[1, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[1, 1].set_ylim([-0.005, 0.2]) 
axes[1, 1].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[1, 1].legend()

plt.tight_layout()
plt.savefig("Imagenes_Resultados/Informacion_Estacionaria_Incoherente_Multiple_I(X_Z)_OR.jpg", dpi=1000)
# %% Graficas multiples de informacion coherentes I(X;Z|Y)

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
fig.suptitle(r"Información mutua en estado estacionario $I(\langle X \rangle;\langle Z \rangle | \langle Y \rangle)$" + "\n" + r"  FFL Incoherentes con diferente coeficiente de Hill |  Compuerta lógica OR", fontsize = 16)

regulacion_K = 0
# Puedes acceder a cada eje individualmente utilizando la notación de índices
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][0]), label = "I1", color = Red_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][1]), label = "I2", color = Orange_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][2]), label = "I3", color = Green_color)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][3]), label = "I4", color = Blue_color)
axes[0, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle | \langle Y \rangle)$", fontsize = 14)
axes[0, 0].set_ylim([-0.005, 0.08]) 
axes[0, 0].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[0, 0].legend() 

regulacion_K = 1
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][0]), label = "I1", color = Red_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][1]), label = "I2", color = Orange_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][2]), label = "I3", color = Green_color)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][3]), label = "I4", color = Blue_color)
axes[0, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle | \langle Y \rangle)$", fontsize = 14)
axes[0, 1].set_ylim([-0.005, 0.08]) 
axes[0, 1].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[0, 1].legend()

regulacion_K = 2
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][0]), label = "I1", color = Red_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][1]), label = "I2", color = Orange_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][2]), label = "I3", color = Green_color)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][3]), label = "I4", color = Blue_color)
axes[1, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle | \langle Y \rangle)$", fontsize = 14)
axes[1, 0].set_ylim([-0.005, 0.08]) 
axes[1, 0].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[1, 0].legend()

regulacion_K = 3
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][0]), label = "I1", color = Red_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][1]), label = "I2", color = Orange_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][2]), label = "I3", color = Green_color)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_Z_X_Y[regulacion_K][3]), label = "I4", color = Blue_color)
axes[1, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle| \langle Y \rangle)$", fontsize = 14)
axes[1, 1].set_ylim([-0.005, 0.08]) 
axes[1, 1].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[1, 1].legend()

plt.tight_layout()

plt.savefig("Imagenes_Resultados/Informacion_Estacionaria_Incoherente_Multiple_I(X_Z|Y)_OR.jpg", dpi=1000)

