#%%
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
from tqdm import tqdm
import math
#%%

def maximizacion_informaciones(simulacion_proteina_X, simulacion_proteina_Y, simulacion_proteina_Z, Posicion_K):
    informacion_maxima_XYZ_cruzados = []
    informacion_maxima_XZ_cruzados = []
    informacion_maxima_YZ_cruzados = []
    informacion_maxima_YZ_mismo_tiempo = []

    tiempo_maximizacion_Y_XYZ = []
    tiempo_maximizacion_Z_XYZ = []
    tiempo_maximizacion_Y_XZ = []
    tiempo_maximizacion_Z_XZ = []

    tiempo_propio_X = np.linspace(9,349,35)
    tiempo_propio_Y = np.linspace(9,349,35)
    tiempo_propio_Z =np.linspace(9,349,35)
    for Tau_X in tqdm(tiempo_propio_X): 
        Tau_X = int(Tau_X)
        informacion_maxima_X_Y_Z = 0
        informacion_maxima_X_Z = 0
        informacion_maxima_Y_Z_mismo_tiempo = 0
        informacion_maxima_Y_Z = 0

        tiempo_maximo_Y_XYZ = Tau_X
        tiempo_maximo_Z_XYZ = Tau_X

        tiempo_maximo_Y_XZ = Tau_X
        tiempo_maximo_Z_XZ = Tau_X

        for Tau_Y in tiempo_propio_Y:
            Tau_Y = int(Tau_Y)
            for Tau_Z in tiempo_propio_Z:
                Tau_Z = int(Tau_Z)
                Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X,simulacion_proteina_Y, simulacion_proteina_Z, Tau_X, Tau_Y, Tau_Z, Posicion_K)
                
                if not  math.isnan(Informacion_Z_X_Y):
                    if Informacion_Z_X_Y>= informacion_maxima_X_Y_Z:
                        informacion_maxima_X_Y_Z = Informacion_Z_X_Y
                        informacion_maxima_Y_Z_mismo_tiempo = Informacion_Z_Y

                        tiempo_maximo_Y_XYZ = Tau_Y
                        tiempo_maximo_Z_XYZ = Tau_Z

                if not  math.isnan(Informacion_Z_X):
                    if Informacion_Z_X>= informacion_maxima_X_Z:
                        informacion_maxima_X_Z = Informacion_Z_X
                        tiempo_maximo_Y_XZ = Tau_Y
                        tiempo_maximo_Z_XZ = Tau_Z

                if not  math.isnan(Informacion_Z_Y):
                    if Informacion_Z_Y>= informacion_maxima_Y_Z:
                        informacion_maxima_Y_Z = Informacion_Z_Y


        informacion_maxima_XYZ_cruzados.append(informacion_maxima_X_Y_Z)
        informacion_maxima_XZ_cruzados.append(informacion_maxima_X_Z)
        informacion_maxima_YZ_cruzados.append(informacion_maxima_Y_Z)
        informacion_maxima_YZ_mismo_tiempo.append(informacion_maxima_Y_Z_mismo_tiempo)


        tiempo_maximizacion_Y_XYZ.append(tiempo_maximo_Y_XYZ)
        tiempo_maximizacion_Z_XYZ.append(tiempo_maximo_Z_XYZ)

        tiempo_maximizacion_Y_XZ.append(tiempo_maximo_Y_XZ)
        tiempo_maximizacion_Z_XZ.append(tiempo_maximo_Z_XZ)

    return informacion_maxima_XYZ_cruzados, informacion_maxima_XZ_cruzados, informacion_maxima_YZ_cruzados, informacion_maxima_YZ_mismo_tiempo, [tiempo_maximizacion_Y_XYZ, tiempo_maximizacion_Z_XYZ], [tiempo_maximizacion_Y_XZ, tiempo_maximizacion_Z_XZ]
#%% Funcion para calcular las informaciones con dependencia temporal
def funcion_informaciones_temporales(simulacion_proteina_X,simulacion_proteina_Y, simulacion_proteina_Z, Tau_X, Tau_Y, Tau_Z, posicion_K):
    data = {'X': simulacion_proteina_X[posicion_K][:, Tau_X].flatten(),
            'Y': simulacion_proteina_Y[posicion_K][:, Tau_Y].flatten(),
            'Z': simulacion_proteina_Z[posicion_K][:, Tau_Z].flatten()}

    df = pd.DataFrame(data)

    cov_matrix = pd.DataFrame.cov(df)
    cov_matrix = np.array(cov_matrix)

    Informacion_Y_X = (1/2)*np.log2((cov_matrix[0][0]* cov_matrix[1][1])/(cov_matrix[0][0]* cov_matrix[1][1] - (cov_matrix[0][1])**2))
    Informacion_Z_X = (1/2)*np.log2((cov_matrix[0][0]* cov_matrix[2][2])/(cov_matrix[0][0]* cov_matrix[2][2] - (cov_matrix[0][2])**2))
    Informacion_Z_Y = (1/2)*np.log2((cov_matrix[1][1]* cov_matrix[2][2])/(cov_matrix[1][1]* cov_matrix[2][2] - (cov_matrix[1][2])**2))
    Informacion_Z_X_Y = (1/2)*np.log2((cov_matrix[0][0]*cov_matrix[1][1]* cov_matrix[2][2] - cov_matrix[0][1]**2)/(np.linalg.det(cov_matrix)))

    return Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y
# %%

#array_nuevo = np.load('Simulacion_FFL_C1_AND.npy', allow_pickle=True).item()
simulacion_FFL_C1 = np.load('/Users/rafaelvelasquez/Documentos/GITHUB/Simulaciones_Resultados/Simulacion_FFL_C1_AND_final.npy', allow_pickle=True).item()
simulacion_FFL_C2 = np.load('/Users/rafaelvelasquez/Documentos/GITHUB/Simulaciones_Resultados/Simulacion_FFL_C2_AND_final.npy', allow_pickle=True).item()
#simulacion_FFL_C3 = np.load('Simulacion_FFL_C3_AND_final.npy', allow_pickle=True).item()
#simulacion_FFL_C4 = np.load('Simulacion_FFL_C4_AND_final.npy', allow_pickle=True).item()

simulacion_FFL_I1 = np.load('/Users/rafaelvelasquez/Documentos/GITHUB/Simulaciones_Resultados/Simulacion_FFL_I1_AND_final.npy', allow_pickle=True).item()
#simulacion_FFL_I2 = np.load('Simulacion_FFL_I2_AND_final.npy', allow_pickle=True).item()
simulacion_FFL_I3 = np.load('/Users/rafaelvelasquez/Documentos/GITHUB/Simulaciones_Resultados/Simulacion_FFL_I3_AND_final.npy', allow_pickle=True).item()
#simulacion_FFL_I4 = np.load('Simulacion_FFL_I4_AND_final.npy', allow_pickle=True).item()
#
#%% SIMULACIONES DIFERENTES FLL

simulacion_proteina_X_C1 = []
for valor_K in range(0,7):
    array_X = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_2"][0][valor_K])
    simulacion_proteina_X_C1.append(array_X)
simulacion_proteina_X_C1 = np.array(simulacion_proteina_X_C1)

simulacion_proteina_Y_C1 = []
for valor_K in range(0,7):
    array_Y = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_2"][1][valor_K])
    simulacion_proteina_Y_C1.append(array_Y)
simulacion_proteina_Y_C1 = np.array(simulacion_proteina_Y_C1)

simulacion_proteina_Z_C1 = []
for valor_K in range(0,7):
    array_Z = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_2"][2][valor_K])
    simulacion_proteina_Z_C1.append(array_Z)
simulacion_proteina_Z_C1 = np.array(simulacion_proteina_Z_C1)

#%%

#%% Simulacion Visualizacion XYZ
Tau_X = 200
Tau_Y = Tau_X
tiempo_propio_Z = np.linspace(0,349,350)
informacion_total = []
for Tau_Y in tqdm(tiempo_propio_Y):
    Tau_Y = int(Tau_Y)
    informacion_X_Y_Z_tiempo_propio_Z = []
    for Tau_Z in tiempo_propio_Z:
        Tau_Z = int(Tau_Z)
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_C1,simulacion_proteina_Y_C1, simulacion_proteina_Z_C1, Tau_X, Tau_Y, Tau_Z, Posicion_K)
        informacion_X_Y_Z_tiempo_propio_Z.append(Informacion_Z_X_Y)
    informacion_total.append(informacion_X_Y_Z_tiempo_propio_Z)
#%%

plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales C1 $\tau_{X}$ fixed value", fontsize = 16)
plt.xlabel(r"Tiempo propio Z", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I (Z;X,Y)$", fontsize = 14)
plt.plot(tiempo_propio_Z, informacion_total[100])




#%%
#____________________________________________________________________________________________________
simulacion_proteina_X_C2 = []
for valor_K in range(0,7):
    array_X = np.array(simulacion_FFL_C2[f"Coeficiente_Hill_2"][0][valor_K])
    simulacion_proteina_X_C2.append(array_X)
simulacion_proteina_X_C2 = np.array(simulacion_proteina_X_C2)

simulacion_proteina_Y_C2 = []
for valor_K in range(0,7):
    array_Y = np.array(simulacion_FFL_C2[f"Coeficiente_Hill_2"][1][valor_K])
    simulacion_proteina_Y_C2.append(array_Y)
simulacion_proteina_Y_C2 = np.array(simulacion_proteina_Y_C2)

simulacion_proteina_Z_C2 = []
for valor_K in range(0,7):
    array_Z = np.array(simulacion_FFL_C2[f"Coeficiente_Hill_2"][2][valor_K])
    simulacion_proteina_Z_C2.append(array_Z)
simulacion_proteina_Z_C2 = np.array(simulacion_proteina_Z_C2)

"""
informacion_Y_X_cruzada_C2 = []
Informacion_Z_X_cruzada_C2 = []
Informacion_Z_Y_cruzada_C2 = []
Informacion_Z_X_Y_cruzada_C2 = []

for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_C2,simulacion_proteina_Y_C2, simulacion_proteina_Z_C2, Tau_Y, Tau_Z, 0)
        
        informacion_Y_X_cruzada_C2.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_C2.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_C2.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_C2.append(Informacion_Z_X_Y)
"""
#%%
tiempos_C2 = []
informacion_maxima_C2 = []
for Tau_X in tqdm(tiempo_propio_X): 
    informacion_Y_X_cruzada_C2 = []
    Informacion_Z_X_cruzada_C2 = []
    Informacion_Z_Y_cruzada_C2 = []
    Informacion_Z_X_Y_cruzada_C2 = []
    for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_C2,simulacion_proteina_Y_C2, simulacion_proteina_Z_C2, Tau_X, Tau_Y, Tau_Z, Posicion_K)
        
        informacion_Y_X_cruzada_C2.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_C2.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_C2.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_C2.append(Informacion_Z_X_Y)

    valor_tiempo_max = np.nanargmax(Informacion_Z_X_cruzada_C2)
    delta_T = valor_tiempo_max - Tau_X

    tiempos_C2.append(delta_T)
    informacion_maxima_C2.append(Informacion_Z_Y_cruzada_C2[valor_tiempo_max])
#%%
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales C2", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z,X)$", fontsize = 14)
plt.plot(tiempo_propio_X[0:300], informacion_maxima_C2[0:300])
#%%
#____________________________________________________________________________________________________
simulacion_proteina_X_C3 = []
for valor_K in range(0,7):
    array_X = np.array(simulacion_FFL_C3[f"Coeficiente_Hill_2"][0][valor_K])
    simulacion_proteina_X_C3.append(array_X)
simulacion_proteina_X_C3 = np.array(simulacion_proteina_X_C3)

simulacion_proteina_Y_C3 = []
for valor_K in range(0,7):
    array_Y = np.array(simulacion_FFL_C3[f"Coeficiente_Hill_2"][1][valor_K])
    simulacion_proteina_Y_C3.append(array_Y)
simulacion_proteina_Y_C3 = np.array(simulacion_proteina_Y_C3)

simulacion_proteina_Z_C3 = []
for valor_K in range(0,7):
    array_Z = np.array(simulacion_FFL_C3[f"Coeficiente_Hill_2"][2][valor_K])
    simulacion_proteina_Z_C3.append(array_Z)
simulacion_proteina_Z_C3 = np.array(simulacion_proteina_Z_C3)
#%%
"""
informacion_Y_X_cruzada_C3 = []
Informacion_Z_X_cruzada_C3 = []
Informacion_Z_Y_cruzada_C3 = []
Informacion_Z_X_Y_cruzada_C3 = []


for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_C3,simulacion_proteina_Y_C3, simulacion_proteina_Z_C3, Tau_Y, Tau_Z, 0)
        
        informacion_Y_X_cruzada_C3.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_C3.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_C3.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_C3.append(Informacion_Z_X_Y)
"""

tiempos_C3 = []
informacion_maxima_C3 = []
for Tau_X in tqdm(tiempo_propio_X): 
    informacion_Y_X_cruzada_C3 = []
    Informacion_Z_X_cruzada_C3 = []
    Informacion_Z_Y_cruzada_C3 = []
    Informacion_Z_X_Y_cruzada_C3 = []
    for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_C3,simulacion_proteina_Y_C3, simulacion_proteina_Z_C3, Tau_X, Tau_Y, Tau_Z, Posicion_K)
        
        informacion_Y_X_cruzada_C3.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_C3.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_C3.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_C3.append(Informacion_Z_X_Y)

    valor_tiempo_max = np.nanargmax(Informacion_Z_X_cruzada_C3)
    delta_T = valor_tiempo_max - Tau_X

    tiempos_C3.append(delta_T)
    informacion_maxima_C3.append(Informacion_Z_Y_cruzada_C3[valor_tiempo_max])
#%%
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales C3", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z,X)$", fontsize = 14)
plt.plot(tiempo_propio_X[0:300], informacion_maxima_C3[0:300])

#%%
#____________________________________________________________________________________________________
simulacion_proteina_X_C4 = []
tiempo_propio_X = np.arange(10,350)
tiempo_propio_Y = np.arange(10,350)
tiempo_propio_Z = np.arange(10,350)
for valor_K in range(0,7):
    array_X = np.array(simulacion_FFL_C4[f"Coeficiente_Hill_2"][0][valor_K])
    simulacion_proteina_X_C4.append(array_X)
simulacion_proteina_X_C4 = np.array(simulacion_proteina_X_C4)

simulacion_proteina_Y_C4 = []
for valor_K in range(0,7):
    array_Y = np.array(simulacion_FFL_C4[f"Coeficiente_Hill_2"][1][valor_K])
    simulacion_proteina_Y_C4.append(array_Y)
simulacion_proteina_Y_C4 = np.array(simulacion_proteina_Y_C4)

simulacion_proteina_Z_C4 = []
for valor_K in range(0,7):
    array_Z = np.array(simulacion_FFL_C4[f"Coeficiente_Hill_2"][2][valor_K])
    simulacion_proteina_Z_C4.append(array_Z)
simulacion_proteina_Z_C4 = np.array(simulacion_proteina_Z_C4)
#%%
"""
informacion_Y_X_cruzada_C4 = []
Informacion_Z_X_cruzada_C4 = []
Informacion_Z_Y_cruzada_C4 = []
Informacion_Z_X_Y_cruzada_C4 = []

for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_C4,simulacion_proteina_Y_C4, simulacion_proteina_Z_C4, Tau_Y, Tau_Z, 0)
        
        informacion_Y_X_cruzada_C4.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_C4.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_C4.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_C4.append(Informacion_Z_X_Y)
"""
tiempos_C4 = []
informacion_maxima_C4 = []
for Tau_X in tqdm(tiempo_propio_X): 
    informacion_Y_X_cruzada_C4 = []
    Informacion_Z_X_cruzada_C4 = []
    Informacion_Z_Y_cruzada_C4 = []
    Informacion_Z_X_Y_cruzada_C4 = []
    for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_C4,simulacion_proteina_Y_C4, simulacion_proteina_Z_C4, Tau_X, Tau_Y, Tau_Z, Posicion_K)
        
        informacion_Y_X_cruzada_C4.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_C4.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_C4.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_C4.append(Informacion_Z_X_Y)

    valor_tiempo_max = np.nanargmax(Informacion_Z_X_Y_cruzada_C4)
    delta_T = valor_tiempo_max - Tau_X

    tiempos_C4.append(delta_T)
    informacion_maxima_C4.append(Informacion_Z_X_Y_cruzada_C4[valor_tiempo_max])

#%% CALCULO DE INFORMACION FINALES
informacion_maxima_XYZ_cruzados_C1, informacion_maxima_XZ_cruzados_C1, informacion_maxima_YZ_cruzados_C1, informacion_maxima_YZ_mismo_tiempo_C1, tiempos_maximizados_XYZ_C1, tiempos_maximizados_XZ_C1 = maximizacion_informaciones(simulacion_proteina_X_C1, simulacion_proteina_Y_C1, simulacion_proteina_Z_C1, Posicion_K = 3)
informacion_maxima_XYZ_cruzados_C2, informacion_maxima_XZ_cruzados_C2, informacion_maxima_YZ_cruzados_C2, informacion_maxima_YZ_mismo_tiempo_C2, tiempos_maximizados_XYZ_C2, tiempos_maximizados_XZ_C2 = maximizacion_informaciones(simulacion_proteina_X_C2, simulacion_proteina_Y_C2, simulacion_proteina_Z_C2, Posicion_K = 3)
#%%
tiempo_propio_X = np.linspace(9,349,35)
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales $I(Z;X,Y)(\tau_x) - I(Z;Y)(\tau_x) $", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z;X,Y)$", fontsize = 14)
plt.plot(tiempo_propio_X, np.array(informacion_maxima_XYZ_cruzados_C1) - np.array(informacion_maxima_YZ_mismo_tiempo_C1), label = r"$I(Z;X,Y)$ C1")
plt.plot(tiempo_propio_X, np.array(informacion_maxima_XYZ_cruzados_C2) - np.array(informacion_maxima_YZ_mismo_tiempo_C2), label = r"$I(Z;X,Y)$ C2")
plt.legend()

tiempo_propio_X = np.linspace(9,349,35)
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales $I(Z;X)(\tau_x)$", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z;X,Y)$", fontsize = 14)

plt.plot(tiempo_propio_X, np.array(informacion_maxima_XZ_cruzados_C1), label = r"$I(Z;X,Y)$ C1")
plt.plot(tiempo_propio_X, np.array(informacion_maxima_XZ_cruzados_C2), label = r"$I(Z;X,Y)$ C2")
#%%



#%%

informacion_maxima_XYZ_cruzados_C1, informacion_maxima_XZ_cruzados_C1, informacion_maxima_YZ_cruzados_C1, informacion_maxima_YZ_mismo_tiempo_C1, tiempos_maximizados_XYZ_C1, tiempos_maximizados_XZ_C1 = maximizacion_informaciones(simulacion_proteina_X_C1, simulacion_proteina_Y_C1, simulacion_proteina_Z_C1, Posicion_K = 3)
informacion_maxima_XYZ_cruzados_C2, informacion_maxima_XZ_cruzados_C2, tiempos_maximizados_XYZ_C2, tiempos_maximizados_XZ_C2 = maximizacion_informaciones(simulacion_proteina_X_C2, simulacion_proteina_Y_C2, simulacion_proteina_Z_C2, Posicion_K = 3)
informacion_maxima_XYZ_cruzados_C3, informacion_maxima_XZ_cruzados_C3, tiempos_maximizados_XYZ_C3, tiempos_maximizados_XZ_C3 = maximizacion_informaciones(simulacion_proteina_X_C3, simulacion_proteina_Y_C3, simulacion_proteina_Z_C3, Posicion_K = 3)
informacion_maxima_XYZ_cruzados_C4, informacion_maxima_XZ_cruzados_C4, tiempos_maximizados_XYZ_C4, tiempos_maximizados_XZ_C4 = maximizacion_informaciones(simulacion_proteina_X_C4, simulacion_proteina_Y_C4, simulacion_proteina_Z_C4, Posicion_K = 3)
#%%
tiempo_propio_X = np.linspace(9,349,341)
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales $I(Z;X,Y)$", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z;X,Y)$", fontsize = 14)
plt.plot(tiempo_propio_X, informacion_maxima_XYZ_cruzados_C1, label = r"$I(Z;X,Y)$ C1")
plt.plot(tiempo_propio_X, informacion_maxima_XYZ_cruzados_C2, label = r"$I(Z;X,Y)$ C2")
plt.plot(tiempo_propio_X, informacion_maxima_XYZ_cruzados_C3, label = r"$I(Z;X,Y)$ C3")
plt.plot(tiempo_propio_X, informacion_maxima_XYZ_cruzados_C4, label = r"$I(Z;X,Y)$ C4")
plt.legend()
#%%
tiempo_propio_X = np.linspace(9,349,341)
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales $I(Z;X)$", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z,X)$", fontsize = 14)
plt.plot(tiempo_propio_X, informacion_maxima_XZ_cruzados_C1, label = r"$I(Z;X)$ C1")
plt.plot(tiempo_propio_X, informacion_maxima_XZ_cruzados_C2, label = r"$I(Z;X)$ C2")
plt.plot(tiempo_propio_X, informacion_maxima_XZ_cruzados_C3, label = r"$I(Z;X)$ C3")
plt.plot(tiempo_propio_X, informacion_maxima_XZ_cruzados_C4, label = r"$I(Z;X)$ C4")
plt.legend()
#%%
#____________________________________________________________________________________________________
tiempo_propio_X = np.arange(10,350)
simulacion_proteina_X_I1 = []
for valor_K in range(0,7):
    array_X = np.array(simulacion_FFL_I1[f"Coeficiente_Hill_2"][0][valor_K])
    simulacion_proteina_X_I1.append(array_X)
simulacion_proteina_X_I1 = np.array(simulacion_proteina_X_I1)

simulacion_proteina_Y_I1 = []
for valor_K in range(0,7):
    array_Y = np.array(simulacion_FFL_I1[f"Coeficiente_Hill_2"][1][valor_K])
    simulacion_proteina_Y_I1.append(array_Y)
simulacion_proteina_Y_I1 = np.array(simulacion_proteina_Y_I1)

simulacion_proteina_Z_I1 = []
for valor_K in range(0,7):
    array_Z = np.array(simulacion_FFL_I1[f"Coeficiente_Hill_2"][2][valor_K])
    simulacion_proteina_Z_I1.append(array_Z)
simulacion_proteina_Z_I1 = np.array(simulacion_proteina_Z_I1)

"""
informacion_Y_X_cruzada_I1 = []
Informacion_Z_X_cruzada_I1 = []
Informacion_Z_Y_cruzada_I1 = []
Informacion_Z_X_Y_cruzada_I1 = []


for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_I1,simulacion_proteina_Y_I1, simulacion_proteina_Z_I1, Tau_Y, Tau_Z, 0)
        
        informacion_Y_X_cruzada_I1.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_I1.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_I1.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_I1.append(Informacion_Z_X_Y)
"""
#%%
tiempos_I1 = []
informacion_maxima_I1 = []
for Tau_X in tqdm(tiempo_propio_X): 
    informacion_Y_X_cruzada_I1 = []
    Informacion_Z_X_cruzada_I1 = []
    Informacion_Z_Y_cruzada_I1 = []
    Informacion_Z_X_Y_cruzada_I1 = []
    for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_I1,simulacion_proteina_Y_I1, simulacion_proteina_Z_I1, Tau_X, Tau_Y, Tau_Z, Posicion_K)
        
        informacion_Y_X_cruzada_I1.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_I1.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_I1.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_I1.append(Informacion_Z_X_Y)

    valor_tiempo_max = np.nanargmax(Informacion_Z_X_Y_cruzada_I1)
    delta_T = valor_tiempo_max - Tau_X

    tiempos_I1.append(delta_T)
    informacion_maxima_I1.append(Informacion_Z_X_Y_cruzada_I1[valor_tiempo_max])
#%%
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales I1", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z,X)$", fontsize = 14)
plt.plot(tiempo_propio_X[0:300], informacion_maxima_I1[0:300])

#%%
#____________________________________________________________________________________________________
simulacion_proteina_X_I2 = []
for valor_K in range(0,7):
    array_X = np.array(simulacion_FFL_I2[f"Coeficiente_Hill_2"][0][valor_K])
    simulacion_proteina_X_I2.append(array_X)
simulacion_proteina_X_I2 = np.array(simulacion_proteina_X_I2)

simulacion_proteina_Y_I2 = []
for valor_K in range(0,7):
    array_Y = np.array(simulacion_FFL_I2[f"Coeficiente_Hill_2"][1][valor_K])
    simulacion_proteina_Y_I2.append(array_Y)
simulacion_proteina_Y_I2 = np.array(simulacion_proteina_Y_I2)

simulacion_proteina_Z_I2 = []
for valor_K in range(0,7):
    array_Z = np.array(simulacion_FFL_I2[f"Coeficiente_Hill_2"][2][valor_K])
    simulacion_proteina_Z_I2.append(array_Z)
simulacion_proteina_Z_I2 = np.array(simulacion_proteina_Z_I2)


"""
informacion_Y_X_cruzada_I2 = []
Informacion_Z_X_cruzada_I2 = []
Informacion_Z_Y_cruzada_I2 = []
Informacion_Z_X_Y_cruzada_I2 = []


for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_I2,simulacion_proteina_Y_I2, simulacion_proteina_Z_I2, Tau_Y, Tau_Z, 0)
        
        informacion_Y_X_cruzada_I2.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_I2.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_I2.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_I2.append(Informacion_Z_X_Y)
"""

tiempos_I2= []
informacion_maxima_I2 = []
for Tau_X in tqdm(tiempo_propio_X): 
    informacion_Y_X_cruzada_I2 = []
    Informacion_Z_X_cruzada_I2 = []
    Informacion_Z_Y_cruzada_I2 = []
    Informacion_Z_X_Y_cruzada_I2 = []
    for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_I2,simulacion_proteina_Y_I2, simulacion_proteina_Z_I2, Tau_X, Tau_Y, Tau_Z, Posicion_K)
        
        informacion_Y_X_cruzada_I2.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_I2.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_I2.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_I2.append(Informacion_Z_X_Y)

    valor_tiempo_max = np.nanargmax(Informacion_Z_X_Y_cruzada_I2)
    delta_T = valor_tiempo_max - Tau_X

    tiempos_I2.append(delta_T)
    informacion_maxima_I2.append(Informacion_Z_X_Y_cruzada_I2[valor_tiempo_max])
#%%
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales I2", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z,X)$", fontsize = 14)
plt.plot(tiempo_propio_X[0:300], informacion_maxima_I2[0:300])


#%%
#____________________________________________________________________________________________________
simulacion_proteina_X_I3 = []
for valor_K in range(0,7):
    array_X = np.array(simulacion_FFL_I3[f"Coeficiente_Hill_2"][0][valor_K])
    simulacion_proteina_X_I3.append(array_X)
simulacion_proteina_X_I3 = np.array(simulacion_proteina_X_I3)

simulacion_proteina_Y_I3 = []
for valor_K in range(0,7):
    array_Y = np.array(simulacion_FFL_I3[f"Coeficiente_Hill_2"][1][valor_K])
    simulacion_proteina_Y_I3.append(array_Y)
simulacion_proteina_Y_I3 = np.array(simulacion_proteina_Y_I3)

simulacion_proteina_Z_I3 = []
for valor_K in range(0,7):
    array_Z = np.array(simulacion_FFL_I3[f"Coeficiente_Hill_2"][2][valor_K])
    simulacion_proteina_Z_I3.append(array_Z)
simulacion_proteina_Z_I3 = np.array(simulacion_proteina_Z_I3)
#%%
"""
informacion_Y_X_cruzada_I3 = []
Informacion_Z_X_cruzada_I3 = []
Informacion_Z_Y_cruzada_I3 = []
Informacion_Z_X_Y_cruzada_I3 = []

for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_I3,simulacion_proteina_Y_I3, simulacion_proteina_Z_I3, Tau_Y, Tau_Z, 0)
        
        informacion_Y_X_cruzada_I3.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_I3.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_I3.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_I3.append(Informacion_Z_X_Y)
"""

tiempos_I3= []
informacion_maxima_I3 = []
for Tau_X in tqdm(tiempo_propio_X): 
    informacion_Y_X_cruzada_I3 = []
    Informacion_Z_X_cruzada_I3 = []
    Informacion_Z_Y_cruzada_I3 = []
    Informacion_Z_X_Y_cruzada_I3 = []
    for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_I3,simulacion_proteina_Y_I3, simulacion_proteina_Z_I3, Tau_X, Tau_Y, Tau_Z, Posicion_K)
        
        informacion_Y_X_cruzada_I3.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_I3.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_I3.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_I3.append(Informacion_Z_X_Y)

    valor_tiempo_max = np.nanargmax(Informacion_Z_X_Y_cruzada_I3)
    delta_T = valor_tiempo_max - Tau_X

    tiempos_I3.append(delta_T)
    informacion_maxima_I3.append(Informacion_Z_X_Y_cruzada_I3[valor_tiempo_max])
#%%
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales I3", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z,X)$", fontsize = 14)
plt.plot(tiempo_propio_X[0:300], informacion_maxima_I3[0:300])

#%%
#____________________________________________________________________________________________________
simulacion_proteina_X_I4 = []
for valor_K in range(0,7):
    array_X = np.array(simulacion_FFL_I4[f"Coeficiente_Hill_2"][0][valor_K])
    simulacion_proteina_X_I4.append(array_X)
simulacion_proteina_X_I4 = np.array(simulacion_proteina_X_I4)

simulacion_proteina_Y_I4 = []
for valor_K in range(0,7):
    array_Y = np.array(simulacion_FFL_I4[f"Coeficiente_Hill_2"][1][valor_K])
    simulacion_proteina_Y_I4.append(array_Y)
simulacion_proteina_Y_I4 = np.array(simulacion_proteina_Y_I4)

simulacion_proteina_Z_I4 = []
for valor_K in range(0,7):
    array_Z = np.array(simulacion_FFL_I4[f"Coeficiente_Hill_2"][2][valor_K])
    simulacion_proteina_Z_I4.append(array_Z)
simulacion_proteina_Z_I4 = np.array(simulacion_proteina_Z_I4)

"""
informacion_Y_X_cruzada_I4 = []
Informacion_Z_X_cruzada_I4 = []
Informacion_Z_Y_cruzada_I4 = []
Informacion_Z_X_Y_cruzada_I4 = []



for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_I4,simulacion_proteina_Y_I4, simulacion_proteina_Z_I4, Tau_Y, Tau_Z, 0)
        
        informacion_Y_X_cruzada_I4.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_I4.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_I4.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_I4.append(Informacion_Z_X_Y)
"""

tiempos_I4= []
informacion_maxima_I4 = []
for Tau_X in tqdm(tiempo_propio_X): 
    informacion_Y_X_cruzada_I4 = []
    Informacion_Z_X_cruzada_I4 = []
    Informacion_Z_Y_cruzada_I4 = []
    Informacion_Z_X_Y_cruzada_I4 = []
    for Tau_Z in tiempo_propio_Z:
        Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X_I4,simulacion_proteina_Y_I4, simulacion_proteina_Z_I4, Tau_X, Tau_Y, Tau_Z, Posicion_K)
        
        informacion_Y_X_cruzada_I4.append(Informacion_Y_X)
        Informacion_Z_X_cruzada_I4.append(Informacion_Z_X)
        Informacion_Z_Y_cruzada_I4.append(Informacion_Z_Y)
        Informacion_Z_X_Y_cruzada_I4.append(Informacion_Z_X_Y)

    valor_tiempo_max = np.nanargmax(Informacion_Z_X_Y_cruzada_I4)
    delta_T = valor_tiempo_max - Tau_X

    tiempos_I4.append(delta_T)
    informacion_maxima_I4.append(Informacion_Z_X_Y_cruzada_I4[valor_tiempo_max])

#%%
informacion_maxima_XYZ_cruzados_I1, informacion_maxima_XZ_cruzados_I1, tiempos_maximizados_XYZ_I1, tiempos_maximizados_XZ_I1 = maximizacion_informaciones(simulacion_proteina_X_I1, simulacion_proteina_Y_I1, simulacion_proteina_Z_I1, Posicion_K = 3)
informacion_maxima_XYZ_cruzados_I2, informacion_maxima_XZ_cruzados_I2, tiempos_maximizados_XYZ_I2, tiempos_maximizados_XZ_I2 = maximizacion_informaciones(simulacion_proteina_X_I2, simulacion_proteina_Y_I2, simulacion_proteina_Z_I2, Posicion_K = 3)
informacion_maxima_XYZ_cruzados_I3, informacion_maxima_XZ_cruzados_I3, tiempos_maximizados_XYZ_I3, tiempos_maximizados_XZ_I3 = maximizacion_informaciones(simulacion_proteina_X_I3, simulacion_proteina_Y_I3, simulacion_proteina_Z_I3, Posicion_K = 3)
informacion_maxima_XYZ_cruzados_I4, informacion_maxima_XZ_cruzados_I4, tiempos_maximizados_XYZ_I4, tiempos_maximizados_XZ_I4 = maximizacion_informaciones(simulacion_proteina_X_I4, simulacion_proteina_Y_I4, simulacion_proteina_Z_I4, Posicion_K = 3)

#%%
informacion_maxima_XYZ_cruzados_I1, informacion_maxima_XZ_cruzados_I1, informacion_maxima_YZ_cruzados_I1, informacion_maxima_YZ_mismo_tiempo_I1, tiempos_maximizados_XYZ_I1, tiempos_maximizados_XZ_I1 = maximizacion_informaciones(simulacion_proteina_X_I1, simulacion_proteina_Y_I1, simulacion_proteina_Z_I1, Posicion_K = 3)
informacion_maxima_XYZ_cruzados_I3, informacion_maxima_XZ_cruzados_I3, informacion_maxima_YZ_cruzados_I3, informacion_maxima_YZ_mismo_tiempo_I3, tiempos_maximizados_XYZ_I3, tiempos_maximizados_XZ_I3 = maximizacion_informaciones(simulacion_proteina_X_I3, simulacion_proteina_Y_I3, simulacion_proteina_Z_I3, Posicion_K = 3)
#%%
tiempo_propio_X = np.linspace(9,349,35)
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales $I(Z;X,Y)(\tau_x) - I(Z;Y)(\tau_x) $", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z;X,Y)$", fontsize = 14)
plt.plot(tiempo_propio_X, np.array(informacion_maxima_XYZ_cruzados_I1) - np.array(informacion_maxima_YZ_mismo_tiempo_I1), label = r"$I(Z;X,Y)$ I1")
plt.plot(tiempo_propio_X, np.array(informacion_maxima_XYZ_cruzados_I3) - np.array(informacion_maxima_YZ_mismo_tiempo_I3), label = r"$I(Z;X,Y)$ I3")
plt.legend()

tiempo_propio_X = np.linspace(9,349,35)
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales $I(Z;X)(\tau_x)$", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z;X,Y)$", fontsize = 14)
plt.plot(tiempo_propio_X, np.array(informacion_maxima_XZ_cruzados_I1), label = r"$I(Z;X,Y)$ I1")
plt.plot(tiempo_propio_X, np.array(informacion_maxima_XZ_cruzados_I3), label = r"$I(Z;X,Y)$ I3")
plt.legend()
#%%
tiempo_propio_X = np.linspace(9,349,35)
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales $I(Z;X,Y)$", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z;X,Y)$", fontsize = 14)
plt.plot(tiempo_propio_X, informacion_maxima_XYZ_cruzados_I1, label = r"$I(Z;X,Y)$ C1")
plt.plot(tiempo_propio_X, informacion_maxima_XYZ_cruzados_I2, label = r"$I(Z;X,Y)$ C2")
plt.plot(tiempo_propio_X, informacion_maxima_XYZ_cruzados_I3, label = r"$I(Z;X,Y)$ C3")
plt.plot(tiempo_propio_X, informacion_maxima_XYZ_cruzados_I4, label = r"$I(Z;X,Y)$ C4")
plt.legend()
#%%
tiempo_propio_X = np.linspace(9,349,35)
plt.figure(figsize=(8,5))
plt.title(r"Informaciones temporales $I(Z;X)$", fontsize = 16)
plt.xlabel(r"Tiempo propio X", fontsize = 14)
plt.ylabel(r"Informacion Mutua $I(Z,X)$", fontsize = 14)
plt.plot(tiempo_propio_X, informacion_maxima_XZ_cruzados_I1, label = r"$I(Z;X)$ C1")
plt.plot(tiempo_propio_X, informacion_maxima_XZ_cruzados_I2, label = r"$I(Z;X)$ C2")
plt.plot(tiempo_propio_X, informacion_maxima_XZ_cruzados_I3, label = r"$I(Z;X)$ C3")
plt.plot(tiempo_propio_X, informacion_maxima_XZ_cruzados_I4, label = r"$I(Z;X)$ C4")
plt.legend()




#%%

maximo_valor_C1 = np.nanargmax(Informacion_Z_X_Y_cruzada_C1)
plt.figure()
plt.title(r"Informaciones temporales C1", fontsize = 16)
plt.xlabel(r"Tiempo propio Z", fontsize = 14)
plt.ylabel(r"Informacion Mutua", fontsize = 14)
plt.plot(tiempo_propio_Z, informacion_Y_X_cruzada_C1, label = "Informacion_Y_X")
plt.plot(tiempo_propio_Z, Informacion_Z_X_cruzada_C1, label = "Informacion_Z_X")
plt.plot(tiempo_propio_Z, Informacion_Z_Y_cruzada_C1, label = "Informacion_Y_Z")
plt.plot(tiempo_propio_Z, Informacion_Z_X_Y_cruzada_C1, label = "Informacion_Z_X_Y")
plt.axvline(x = Tau_Y,  color = 'b', label = 'Tiempo propio Y', linestyle = "--")
plt.axvline(x = tiempo_propio_Z[maximo_valor_C1], ymin = 0,  ymax = Informacion_Z_X_Y_cruzada_C1[maximo_valor_C1] ,  color = 'violet', label = 'Valor_maximo', linestyle = "--")
plt.legend()
print(maximo_valor_C1 - Tau_Y)


maximo_valor = np.nanargmax(Informacion_Z_X_Y_cruzada_C2)
plt.figure()
plt.title(r"Informaciones temporales C2", fontsize = 16)
plt.xlabel(r"Tiempo propio Z", fontsize = 14)
plt.ylabel(r"Informacion Mutua", fontsize = 14)
plt.plot(tiempo_propio_Z, informacion_Y_X_cruzada_C2, label = "Informacion_Y_X")
plt.plot(tiempo_propio_Z, Informacion_Z_X_cruzada_C2, label = "Informacion_Z_X")
plt.plot(tiempo_propio_Z, Informacion_Z_Y_cruzada_C2, label = "Informacion_Y_Z")
plt.plot(tiempo_propio_Z, Informacion_Z_X_Y_cruzada_C2, label = "Informacion_Z_X_Y")
plt.axvline(x = Tau_Y, ymin = 0,  ymax = Informacion_Z_X_Y_cruzada_C2[Tau_Y] ,  color = 'b', label = 'Tiempo propio Y', linestyle = "--")
plt.axvline(x = tiempo_propio_Z[maximo_valor], ymin = 0,  ymax = Informacion_Z_X_Y_cruzada_C2[maximo_valor] ,  color = 'violet', label = 'Valor_maximo', linestyle = "--")
plt.legend()
print(maximo_valor - Tau_Y)


maximo_valor = np.nanargmax(Informacion_Z_X_Y_cruzada_C3)
plt.figure()
plt.title(r"Informaciones temporales C3", fontsize = 16)
plt.xlabel(r"Tiempo propio Z", fontsize = 14)
plt.ylabel(r"Informacion Mutua", fontsize = 14)
plt.plot(tiempo_propio_Z, informacion_Y_X_cruzada_C3, label = "Informacion_Y_X")
plt.plot(tiempo_propio_Z, Informacion_Z_X_cruzada_C3, label = "Informacion_Z_X")
plt.plot(tiempo_propio_Z, Informacion_Z_Y_cruzada_C3, label = "Informacion_Y_Z")
plt.plot(tiempo_propio_Z, Informacion_Z_X_Y_cruzada_C3, label = "Informacion_Z_X_Y")
plt.axvline(x = Tau_Y ,  color = 'b', label = 'Tiempo propio Y', linestyle = "--")
plt.axvline(x = tiempo_propio_Z[maximo_valor] ,  color = 'violet', label = 'Valor_maximo', linestyle = "--")
plt.legend()
print(maximo_valor - Tau_Y)


maximo_valor = np.nanargmax(Informacion_Z_X_Y_cruzada_C4)
plt.figure()
plt.title(r"Informaciones temporales C4", fontsize = 16)
plt.xlabel(r"Tiempo propio Z", fontsize = 14)
plt.ylabel(r"Informacion Mutua", fontsize = 14)
plt.plot(tiempo_propio_Z, informacion_Y_X_cruzada_C4, label = "Informacion_Y_X")
plt.plot(tiempo_propio_Z, Informacion_Z_X_cruzada_C4, label = "Informacion_Z_X")
plt.plot(tiempo_propio_Z, Informacion_Z_Y_cruzada_C4, label = "Informacion_Y_Z")
plt.plot(tiempo_propio_Z, Informacion_Z_X_Y_cruzada_C4, label = "Informacion_Z_X_Y")
plt.axvline(x = Tau_Y ,  color = 'b', label = 'Tiempo propio Y', linestyle = "--")
plt.axvline(x = tiempo_propio_Z[maximo_valor],  color = 'violet', label = 'Valor_maximo', linestyle = "--")
plt.legend()
print(maximo_valor - Tau_Y)


maximo_valor = np.nanargmax(Informacion_Z_X_Y_cruzada_I1)
plt.figure()
plt.title(r"Informaciones temporales I1", fontsize = 16)
plt.xlabel(r"Tiempo propio Z", fontsize = 14)
plt.ylabel(r"Informacion Mutua", fontsize = 14)
plt.plot(tiempo_propio_Z, informacion_Y_X_cruzada_I1, label = "Informacion_Y_X")
plt.plot(tiempo_propio_Z, Informacion_Z_X_cruzada_I1, label = "Informacion_Z_X")
plt.plot(tiempo_propio_Z, Informacion_Z_Y_cruzada_I1, label = "Informacion_Y_Z")
plt.plot(tiempo_propio_Z, Informacion_Z_X_Y_cruzada_I1, label = "Informacion_Z_X_Y")
plt.axvline(x = Tau_Y,  color = 'b', label = 'Tiempo propio Y', linestyle = "--")
plt.axvline(x = tiempo_propio_Z[maximo_valor] ,  color = 'violet', label = 'Valor_maximo', linestyle = "--")
plt.legend()
print(maximo_valor - Tau_Y)


maximo_valor = np.nanargmax(Informacion_Z_X_Y_cruzada_I2)
plt.figure()
plt.title(r"Informaciones temporales I2", fontsize = 16)
plt.xlabel(r"Tiempo propio Z", fontsize = 14)
plt.ylabel(r"Informacion Mutua", fontsize = 14)
plt.plot(tiempo_propio_Z, informacion_Y_X_cruzada_I2, label = "Informacion_Y_X")
plt.plot(tiempo_propio_Z, Informacion_Z_X_cruzada_I2, label = "Informacion_Z_X")
plt.plot(tiempo_propio_Z, Informacion_Z_Y_cruzada_I2, label = "Informacion_Y_Z")
plt.plot(tiempo_propio_Z, Informacion_Z_X_Y_cruzada_I2, label = "Informacion_Z_X_Y")
plt.axvline(x = Tau_Y,  color = 'b', label = 'Tiempo propio Y', linestyle = "--")
plt.axvline(x = tiempo_propio_Z[maximo_valor] ,  color = 'violet', label = 'Valor_maximo', linestyle = "--")
plt.legend()

print(maximo_valor - Tau_Y)


maximo_valor = np.nanargmax(Informacion_Z_X_Y_cruzada_I3)
plt.figure()
plt.title(r"Informaciones temporales I3", fontsize = 16)
plt.xlabel(r"Tiempo propio Z", fontsize = 14)
plt.ylabel(r"Informacion Mutua", fontsize = 14)
plt.plot(tiempo_propio_Z, informacion_Y_X_cruzada_I3, label = "Informacion_Y_X")
plt.plot(tiempo_propio_Z, Informacion_Z_X_cruzada_I3, label = "Informacion_Z_X")
plt.plot(tiempo_propio_Z, Informacion_Z_Y_cruzada_I3, label = "Informacion_Y_Z")
plt.plot(tiempo_propio_Z, Informacion_Z_X_Y_cruzada_I3, label = "Informacion_Z_X_Y")
plt.axvline(x = Tau_Y,  color = 'b', label = 'Tiempo propio Y', linestyle = "--")
plt.axvline(x = tiempo_propio_Z[maximo_valor] ,  color = 'violet', label = 'Valor_maximo', linestyle = "--")
plt.legend()

print(maximo_valor - Tau_Y)


maximo_valor = np.nanargmax(Informacion_Z_X_Y_cruzada_I4)
plt.figure()
plt.title(r"Informaciones temporales I4", fontsize = 16)
plt.xlabel(r"Tiempo propio Z", fontsize = 14)
plt.ylabel(r"Informacion Mutua", fontsize = 14)
plt.plot(tiempo_propio_Z, informacion_Y_X_cruzada_I4, label = "Informacion_Y_X")
plt.plot(tiempo_propio_Z, Informacion_Z_X_cruzada_I4, label = "Informacion_Z_X")
plt.plot(tiempo_propio_Z, Informacion_Z_Y_cruzada_I4, label = "Informacion_Y_Z")
plt.plot(tiempo_propio_Z, Informacion_Z_X_Y_cruzada_I4, label = "Informacion_Z_X_Y")
plt.axvline(x = Tau_Y,  color = 'b', label = 'Tiempo propio Y', linestyle = "--")
plt.axvline(x = tiempo_propio_Z[maximo_valor] ,  color = 'violet', label = 'Valor_maximo', linestyle = "--")
plt.legend()

print(maximo_valor - Tau_Y)

# %%
