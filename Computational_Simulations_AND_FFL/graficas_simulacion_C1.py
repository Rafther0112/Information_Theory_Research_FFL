#%%
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
# %%
#array_nuevo = np.load('Simulacion_FFL_C1_AND.npy', allow_pickle=True).item()
simulacion_FFL_C1 = np.load('Simulacion_FFL_C1_AND_final.npy', allow_pickle=True).item()

def funcion_informaciones_temporales(simulacion_proteina_X,simulacion_proteina_Y, simulacion_proteina_Z, Tau_X, Tau_Y, Tau_Z, posicion_K):
    data = {'X': simulacion_proteina_X[posicion_K][:, Tau_X],
            'Y': simulacion_proteina_Y[posicion_K][:, Tau_Y],
            'Z': simulacion_proteina_Z[posicion_K][:, Tau_Z]}

    df = pd.DataFrame(data)

    cov_matrix = pd.DataFrame.cov(df)
    cov_matrix = np.array(cov_matrix)

    Informacion_Y_X = (1/2)*np.log2((cov_matrix[0][0]* cov_matrix[1][1])/(cov_matrix[0][0]* cov_matrix[1][1] - (cov_matrix[0][1])**2))
    Informacion_Z_X = (1/2)*np.log2((cov_matrix[0][0]* cov_matrix[2][2])/(cov_matrix[0][0]* cov_matrix[2][2] - (cov_matrix[0][2])**2))
    Informacion_Z_Y = (1/2)*np.log2((cov_matrix[1][1]* cov_matrix[2][2])/(cov_matrix[1][1]* cov_matrix[2][2] - (cov_matrix[1][2])**2))
    Informacion_Z_X_Y = (1/2)*np.log2((cov_matrix[0][0]*cov_matrix[1][1]* cov_matrix[2][2] - cov_matrix[0][1]**2)/(np.linalg.det(cov_matrix)))

    return Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y

#%%
simulacion_proteina_X = []
for valor_K in range(0,7):
    array_X = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_2"][0][valor_K])
    simulacion_proteina_X.append(array_X)
simulacion_proteina_X = np.array(simulacion_proteina_X)

simulacion_proteina_Y = []
for valor_K in range(0,7):
    array_Y = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_2"][1][valor_K])
    simulacion_proteina_Y.append(array_Y)
simulacion_proteina_Y = np.array(simulacion_proteina_Y)

simulacion_proteina_Z = []
for valor_K in range(0,7):
    array_Z = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_2"][2][valor_K])
    simulacion_proteina_Z.append(array_Z)
simulacion_proteina_Z = np.array(simulacion_proteina_Z)
#%%
tiempo_propio_X =  np.arange(0,350)
tiempo_propio_Y = np.arange(0,350)
tiempo_propio_Z = np.arange(0,350)
#%%

tiempo_Z_maximizacion = []
informacion_X_Z_maximizada = []
delay_de_maximizacion = []

from tqdm import tqdm 
for Tau_X  in tqdm(tiempo_propio_X):
    informacion_XZ_propia = []
    Tau_Y = -1
    posicion_K = 0
    for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X,simulacion_proteina_Y, simulacion_proteina_Z, Tau_X, Tau_Y, Tau_Z, posicion_K)
        informacion_XZ_propia.append(Informacion_Z_X)
    if not (np.isnan(informacion_XZ_propia).all()):
        tiempo_maximo_Z = np.nanargmax(informacion_XZ_propia)
        informacion_X_Z_maximizada.append(informacion_XZ_propia[tiempo_maximo_Z])
        tiempo_Z_maximizacion.append(tiempo_maximo_Z)
        delay = tiempo_maximo_Z - Tau_X
        delay_de_maximizacion.append(delay)

    else: 
        tiempo_maximo_Z = 0

#%%
plt.figure()
plt.title(r"Informaciones temporales C1", fontsize = 16)
plt.xlabel(r"Tiempo propio Z", fontsize = 14)
plt.ylabel(r"Informacion Mutua", fontsize = 14)
plt.scatter(tiempo_Z_maximizacion, informacion_X_Z_maximizada, label = "Informacion_Z_X")
plt.legend()
#%%
plt.plot(delay_de_maximizacion)
# %%
delay_de_maximizacion
# %%
