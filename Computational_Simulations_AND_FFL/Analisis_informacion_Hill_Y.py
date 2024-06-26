"""
Este Script es para hacer el analisis del experimento de modificar
el valor de la constante de regulación del gen Y sobre el gen Z y 
mirar que variaciones hay sobre la informacion I(Z;X)
"""

#%% HACEMOS LOS IMPORTS DE LAS LIBRERIAS QUE VAMOS A UTILIZAR
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
from tqdm import tqdm
import math
#%% HACEMOS LA LECTURA DE LAS SIMULACIONES

#array_nuevo = np.load('Simulacion_FFL_C1_AND.npy', allow_pickle=True).item()
simulacion_FFL_C1 = np.load('Simulacion_FFL_C1_AND_Modificacion_Hill_Y_inicial.npy', allow_pickle=True).item()
simulacion_FFL_C2 = np.load('Simulacion_FFL_C2_AND_Modificacion_Hill_Y_inicial.npy', allow_pickle=True).item()

simulacion_FFL_I1 = np.load('Simulacion_FFL_I1_AND_Modificacion_Hill_Y_inicial.npy', allow_pickle=True).item()
simulacion_FFL_I3 = np.load('Simulacion_FFL_I3_AND_Modificacion_Hill_Y_inicial.npy', allow_pickle=True).item()

# %%
simulacion_proteina_X_C1 = []
for valor_K in range(0,5):
    array_X = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_3"][0][valor_K])
    simulacion_proteina_X_C1.append(array_X)
simulacion_proteina_X_C1 = np.array(simulacion_proteina_X_C1)

simulacion_proteina_Y_C1 = []
for valor_K in range(0,5):
    array_Y = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_3"][1][valor_K])
    simulacion_proteina_Y_C1.append(array_Y)
simulacion_proteina_Y_C1 = np.array(simulacion_proteina_Y_C1)

simulacion_proteina_Z_C1 = []
for valor_K in range(0,5):
    array_Z = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_3"][2][valor_K])
    simulacion_proteina_Z_C1.append(array_Z)
simulacion_proteina_Z_C1 = np.array(simulacion_proteina_Z_C1)


simulacion_proteina_X_C2 = []
for valor_K in range(0,5):
    array_X = np.array(simulacion_FFL_C2[f"Coeficiente_Hill_3"][0][valor_K])
    simulacion_proteina_X_C2.append(array_X)
simulacion_proteina_X_C2 = np.array(simulacion_proteina_X_C2)

simulacion_proteina_Y_C2 = []
for valor_K in range(0,5):
    array_Y = np.array(simulacion_FFL_C2[f"Coeficiente_Hill_3"][1][valor_K])
    simulacion_proteina_Y_C2.append(array_Y)
simulacion_proteina_Y_C2 = np.array(simulacion_proteina_Y_C2)

simulacion_proteina_Z_C2 = []
for valor_K in range(0,5):
    array_Z = np.array(simulacion_FFL_C2[f"Coeficiente_Hill_3"][2][valor_K])
    simulacion_proteina_Z_C2.append(array_Z)
simulacion_proteina_Z_C2 = np.array(simulacion_proteina_Z_C2)


simulacion_proteina_X_I1 = []
for valor_K in range(0,5):
    array_X = np.array(simulacion_FFL_I1[f"Coeficiente_Hill_3"][0][valor_K])
    simulacion_proteina_X_I1.append(array_X)
simulacion_proteina_X_I1 = np.array(simulacion_proteina_X_I1)

simulacion_proteina_Y_I1 = []
for valor_K in range(0,5):
    array_Y = np.array(simulacion_FFL_I1[f"Coeficiente_Hill_3"][1][valor_K])
    simulacion_proteina_Y_I1.append(array_Y)
simulacion_proteina_Y_I1 = np.array(simulacion_proteina_Y_I1)

simulacion_proteina_Z_I1 = []
for valor_K in range(0,5):
    array_Z = np.array(simulacion_FFL_I1[f"Coeficiente_Hill_3"][2][valor_K])
    simulacion_proteina_Z_I1.append(array_Z)
simulacion_proteina_Z_I1 = np.array(simulacion_proteina_Z_I1)


simulacion_proteina_X_I3 = []
for valor_K in range(0,5):
    array_X = np.array(simulacion_FFL_I3[f"Coeficiente_Hill_3"][0][valor_K])
    simulacion_proteina_X_I3.append(array_X)
simulacion_proteina_X_I3 = np.array(simulacion_proteina_X_I3)

simulacion_proteina_Y_I3 = []
for valor_K in range(0,5):
    array_Y = np.array(simulacion_FFL_I3[f"Coeficiente_Hill_3"][1][valor_K])
    simulacion_proteina_Y_I3.append(array_Y)
simulacion_proteina_Y_I3 = np.array(simulacion_proteina_Y_I3)

simulacion_proteina_Z_I3 = []
for valor_K in range(0,5):
    array_Z = np.array(simulacion_FFL_I3[f"Coeficiente_Hill_3"][2][valor_K])
    simulacion_proteina_Z_I3.append(array_Z)
simulacion_proteina_Z_I3 = np.array(simulacion_proteina_Z_I3)

# %%
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

def maximizacion_informaciones(simulacion_proteina_X, simulacion_proteina_Y, simulacion_proteina_Z, Posicion_K):

    tiempo_propio_X = np.linspace(9,349,35)
    tiempo_propio_Y = np.linspace(9,349,341)
    tiempo_propio_Z =np.linspace(9,349,341)
    matrix_informacion_YX = np.empty(shape=(len(tiempo_propio_X), len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='object')
    matrix_informacion_ZX = np.empty(shape=(len(tiempo_propio_X), len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='object')
    matrix_informacion_ZY = np.empty(shape=(len(tiempo_propio_X), len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='object')
    matrix_informacion_ZXY = np.empty(shape=(len(tiempo_propio_X), len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='object')
    for posicion_X, Tau_X in enumerate(tqdm(tiempo_propio_X)): 
        Tau_X = int(Tau_X)
        for posicion_Y, Tau_Y in enumerate(tiempo_propio_Y):
            Tau_Y = int(Tau_Y)
            for poscion_Z , Tau_Z in enumerate(tiempo_propio_Z):
                Tau_Z = int(Tau_Z)
                Informacion_Y_X, Informacion_Z_X, Informacion_Z_Y, Informacion_Z_X_Y = funcion_informaciones_temporales(simulacion_proteina_X,simulacion_proteina_Y, simulacion_proteina_Z, Tau_X, Tau_Y, Tau_Z, Posicion_K)
                matrix_informacion_YX[posicion_X][posicion_Y][poscion_Z] = Informacion_Y_X
                matrix_informacion_ZX[posicion_X][posicion_Y][poscion_Z] = Informacion_Z_X
                matrix_informacion_ZY[posicion_X][posicion_Y][poscion_Z] = Informacion_Z_Y
                matrix_informacion_ZXY[posicion_X][posicion_Y][poscion_Z] = Informacion_Z_X_Y
                
    return matrix_informacion_YX, matrix_informacion_ZX, matrix_informacion_ZY, matrix_informacion_ZXY
# %%
matrix_info_YX_C1,matrix_info_ZX_C1,matrix_info_ZY_C1,matrix_info_ZXY_C1 = maximizacion_informaciones(simulacion_proteina_X_C1, simulacion_proteina_Y_C1, simulacion_proteina_Z_C1, Posicion_K= 4)
matrix_info_YX_C2,matrix_info_ZX_C2,matrix_info_ZY_C2,matrix_info_ZXY_C2 = maximizacion_informaciones(simulacion_proteina_X_C2, simulacion_proteina_Y_C2, simulacion_proteina_Z_C2, Posicion_K= 4)

matrix_info_YX_I1,matrix_info_ZX_I1,matrix_info_ZY_I1,matrix_info_ZXY_I1 = maximizacion_informaciones(simulacion_proteina_X_I1, simulacion_proteina_Y_I1, simulacion_proteina_Z_I1, Posicion_K= 4)
matrix_info_YX_I3,matrix_info_ZX_I3,matrix_info_ZY_I3,matrix_info_ZXY_I3 = maximizacion_informaciones(simulacion_proteina_X_I3, simulacion_proteina_Y_I3, simulacion_proteina_Z_I3, Posicion_K= 4)

# %%
informacion_total_total_C1 = []
informacion_total_total_C2 = []
#informacion_total_total_C3 = []
#informacion_total_total_C4 = []

informacion_total_total_I1 = []
#informacion_total_total_I2 = []
informacion_total_total_I3 = []
#informacion_total_total_I4 = []

for tiempo_prop_x in tqdm(range(0,35)):

# Assuming you have four matrices of shape (35, 35)
    matrix_dataC1 = matrix_info_ZX_C1[tiempo_prop_x].astype(float)
    matrix_dataC2 = matrix_info_ZX_C2[tiempo_prop_x].astype(float)
    #matrix_dataC3 = matrix_info_ZXY_C3[tiempo_prop_x].astype(float)
    #matrix_dataC4 = matrix_info_ZXY_C4[tiempo_prop_x].astype(float)

    matrix_dataI1 = matrix_info_ZX_I1[tiempo_prop_x].astype(float)
    #matrix_dataI2 = matrix_info_ZXY_I2[tiempo_prop_x].astype(float)
    matrix_dataI3 = matrix_info_ZX_I3[tiempo_prop_x].astype(float)
    #matrix_dataI4 = matrix_info_ZXY_I4[tiempo_prop_x].astype(float)

    average_intensity_C1 = np.nanmean(matrix_dataC1)
    average_intensity_C2 = np.nanmean(matrix_dataC2)
    #average_intensity_C3 = np.mean(matrix_dataC3)
    #average_intensity_C4 = np.mean(matrix_dataC4)

    average_intensity_I1 = np.nanmean(matrix_dataI1)
    #average_intensity_I2 = np.mean(matrix_dataI2)
    average_intensity_I3 = np.nanmean(matrix_dataI3)
    #average_intensity_I4 = np.mean(matrix_dataI4)

    informacion_total_total_C1.append(average_intensity_C1)
    informacion_total_total_C2.append(average_intensity_C2)
    #informacion_total_total_C3.append(average_intensity_C3)
    #informacion_total_total_C4.append(average_intensity_C4)

    informacion_total_total_I1.append(average_intensity_I1)
    #informacion_total_total_I2.append(average_intensity_I2)
    informacion_total_total_I3.append(average_intensity_I3)
    #informacion_total_total_I4.append(average_intensity_I4)
#%%
plt.title(r"Información promedio $I(Z;X)$ FFL C n = 1 | K = 2")
plt.xlabel(r"Tiempo propio X $\tau_x$")
plt.ylabel(r"Informacóon promedio $\langle I(Z;X) \rangle$")

plt.plot(range(0,35), informacion_total_total_C1, label = "FFL C1")
plt.plot(range(0,35), informacion_total_total_C2, label = "FFL C2")
#plt.plot(range(0,70), informacion_total_total_C3, label = "FFL C3")
#plt.plot(range(0,70), informacion_total_total_C4, label = "FFL C4")
plt.legend()
# %%
plt.title(r"Información promedio $I(Z;X)$ FFL I n = 1 | K = 3")
plt.xlabel(r"Tiempo propio X $\tau_x$")
plt.ylabel(r"Informacóon promedio $\langle I(Z;X) \rangle$")
plt.plot(range(0,35), informacion_total_total_I1, label = "FFL I1")
#plt.plot(range(0,70), informacion_total_total_I2, label = "FFL I2")
plt.plot(range(0,35), informacion_total_total_I3, label = "FFL I3")
#plt.plot(range(0,70), informacion_total_total_I4, label = "FFL I4")
plt.legend()
# %%


diccionario_global_FFL_C1 = {}
diccionario_global_FFL_C2 = {}
diccionario_global_FFL_I1 = {}
diccionario_global_FFL_I3 = {}

valores_posibles_Hill = [1,2,3,4]
valores_posibles_Kx = [0,1,2,3,4]
for Hill_modulo_Y in valores_posibles_Hill:

    matrix_info_YX_C1_GLOBAL,matrix_info_ZX_C1_GLOBAL,matrix_info_ZY_C1_GLOBAL,matrix_info_ZXY_C1_GLOBAL = [],[],[],[]
    matrix_info_YX_C2_GLOBAL,matrix_info_ZX_C2_GLOBAL,matrix_info_ZY_C2_GLOBAL,matrix_info_ZXY_C2_GLOBAL = [],[],[],[]

    matrix_info_YX_I1_GLOBAL,matrix_info_ZX_I1_GLOBAL,matrix_info_ZY_I1_GLOBAL,matrix_info_ZXY_I1_GLOBAL = [],[],[],[]
    matrix_info_YX_I3_GLOBAL,matrix_info_ZX_I3_GLOBAL,matrix_info_ZY_I3_GLOBAL,matrix_info_ZXY_I3_GLOBAL = [],[],[],[]

    simulacion_proteina_X_C1 = []
    for valor_K in range(0,5):
        array_X = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_{Hill_modulo_Y}"][0][valor_K])
        simulacion_proteina_X_C1.append(array_X)
    simulacion_proteina_X_C1 = np.array(simulacion_proteina_X_C1)

    simulacion_proteina_Y_C1 = []
    for valor_K in range(0,5):
        array_Y = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_{Hill_modulo_Y}"][1][valor_K])
        simulacion_proteina_Y_C1.append(array_Y)
    simulacion_proteina_Y_C1 = np.array(simulacion_proteina_Y_C1)

    simulacion_proteina_Z_C1 = []
    for valor_K in range(0,5):
        array_Z = np.array(simulacion_FFL_C1[f"Coeficiente_Hill_{Hill_modulo_Y}"][2][valor_K])
        simulacion_proteina_Z_C1.append(array_Z)
    simulacion_proteina_Z_C1 = np.array(simulacion_proteina_Z_C1)


    simulacion_proteina_X_C2 = []
    for valor_K in range(0,5):
        array_X = np.array(simulacion_FFL_C2[f"Coeficiente_Hill_{Hill_modulo_Y}"][0][valor_K])
        simulacion_proteina_X_C2.append(array_X)
    simulacion_proteina_X_C2 = np.array(simulacion_proteina_X_C2)

    simulacion_proteina_Y_C2 = []
    for valor_K in range(0,5):
        array_Y = np.array(simulacion_FFL_C2[f"Coeficiente_Hill_{Hill_modulo_Y}"][1][valor_K])
        simulacion_proteina_Y_C2.append(array_Y)
    simulacion_proteina_Y_C2 = np.array(simulacion_proteina_Y_C2)

    simulacion_proteina_Z_C2 = []
    for valor_K in range(0,5):
        array_Z = np.array(simulacion_FFL_C2[f"Coeficiente_Hill_{Hill_modulo_Y}"][2][valor_K])
        simulacion_proteina_Z_C2.append(array_Z)
    simulacion_proteina_Z_C2 = np.array(simulacion_proteina_Z_C2)


    simulacion_proteina_X_I1 = []
    for valor_K in range(0,5):
        array_X = np.array(simulacion_FFL_I1[f"Coeficiente_Hill_{Hill_modulo_Y}"][0][valor_K])
        simulacion_proteina_X_I1.append(array_X)
    simulacion_proteina_X_I1 = np.array(simulacion_proteina_X_I1)

    simulacion_proteina_Y_I1 = []
    for valor_K in range(0,5):
        array_Y = np.array(simulacion_FFL_I1[f"Coeficiente_Hill_{Hill_modulo_Y}"][1][valor_K])
        simulacion_proteina_Y_I1.append(array_Y)
    simulacion_proteina_Y_I1 = np.array(simulacion_proteina_Y_I1)

    simulacion_proteina_Z_I1 = []
    for valor_K in range(0,5):
        array_Z = np.array(simulacion_FFL_I1[f"Coeficiente_Hill_{Hill_modulo_Y}"][2][valor_K])
        simulacion_proteina_Z_I1.append(array_Z)
    simulacion_proteina_Z_I1 = np.array(simulacion_proteina_Z_I1)


    simulacion_proteina_X_I3 = []
    for valor_K in range(0,5):
        array_X = np.array(simulacion_FFL_I3[f"Coeficiente_Hill_{Hill_modulo_Y}"][0][valor_K])
        simulacion_proteina_X_I3.append(array_X)
    simulacion_proteina_X_I3 = np.array(simulacion_proteina_X_I3)

    simulacion_proteina_Y_I3 = []
    for valor_K in range(0,5):
        array_Y = np.array(simulacion_FFL_I3[f"Coeficiente_Hill_{Hill_modulo_Y}"][1][valor_K])
        simulacion_proteina_Y_I3.append(array_Y)
    simulacion_proteina_Y_I3 = np.array(simulacion_proteina_Y_I3)

    simulacion_proteina_Z_I3 = []
    for valor_K in range(0,5):
        array_Z = np.array(simulacion_FFL_I3[f"Coeficiente_Hill_{Hill_modulo_Y}"][2][valor_K])
        simulacion_proteina_Z_I3.append(array_Z)
    simulacion_proteina_Z_I3 = np.array(simulacion_proteina_Z_I3)


    for Kx in tqdm(valores_posibles_Kx):
        matrix_info_YX_C1,matrix_info_ZX_C1,matrix_info_ZY_C1,matrix_info_ZXY_C1 = maximizacion_informaciones(simulacion_proteina_X_C1, simulacion_proteina_Y_C1, simulacion_proteina_Z_C1, Posicion_K= Kx)
        matrix_info_YX_C2,matrix_info_ZX_C2,matrix_info_ZY_C2,matrix_info_ZXY_C2 = maximizacion_informaciones(simulacion_proteina_X_C2, simulacion_proteina_Y_C2, simulacion_proteina_Z_C2, Posicion_K= Kx)

        matrix_info_YX_I1,matrix_info_ZX_I1,matrix_info_ZY_I1,matrix_info_ZXY_I1 = maximizacion_informaciones(simulacion_proteina_X_I1, simulacion_proteina_Y_I1, simulacion_proteina_Z_I1, Posicion_K= Kx)
        matrix_info_YX_I3,matrix_info_ZX_I3,matrix_info_ZY_I3,matrix_info_ZXY_I3 = maximizacion_informaciones(simulacion_proteina_X_I3, simulacion_proteina_Y_I3, simulacion_proteina_Z_I3, Posicion_K= Kx)

        matrix_info_YX_C1_GLOBAL.append(matrix_info_YX_C1)
        matrix_info_ZX_C1_GLOBAL.append(matrix_info_ZX_C1)
        matrix_info_ZY_C1_GLOBAL.append(matrix_info_ZY_C1)
        matrix_info_ZXY_C1_GLOBAL.append(matrix_info_ZXY_C1)

        matrix_info_YX_C2_GLOBAL.append(matrix_info_YX_C2)
        matrix_info_ZX_C2_GLOBAL.append(matrix_info_ZX_C2)
        matrix_info_ZY_C2_GLOBAL.append(matrix_info_ZY_C2)
        matrix_info_ZXY_C2_GLOBAL.append(matrix_info_ZXY_C2)

        matrix_info_YX_I1_GLOBAL.append(matrix_info_YX_I1)
        matrix_info_ZX_I1_GLOBAL.append(matrix_info_ZX_I1)
        matrix_info_ZY_I1_GLOBAL.append(matrix_info_ZY_I1)
        matrix_info_ZXY_I1_GLOBAL.append(matrix_info_ZXY_I1)

        matrix_info_YX_I3_GLOBAL.append(matrix_info_YX_I3)
        matrix_info_ZX_I3_GLOBAL.append(matrix_info_ZX_I3)
        matrix_info_ZY_I3_GLOBAL.append(matrix_info_ZY_I3)
        matrix_info_ZXY_I3_GLOBAL.append(matrix_info_ZXY_I3)

    diccionario_global_FFL_C1[f"Coeficiente_Hill_{Hill_modulo_Y}"] = [matrix_info_YX_C1_GLOBAL,matrix_info_ZX_C1_GLOBAL,matrix_info_ZY_C1_GLOBAL,matrix_info_ZXY_C1_GLOBAL]
    np.save('MATRICES_INFORMACION_FFL_C1_AND_Modificacion_Hill_Y_inicial.npy', diccionario_global_FFL_C1)
    diccionario_global_FFL_C2[f"Coeficiente_Hill_{Hill_modulo_Y}"] = [matrix_info_YX_C2_GLOBAL,matrix_info_ZX_C2_GLOBAL,matrix_info_ZY_C2_GLOBAL,matrix_info_ZXY_C2_GLOBAL]
    np.save('MATRICES_INFORMACION_FFL_C2_AND_Modificacion_Hill_Y_inicial.npy', diccionario_global_FFL_C2)
    
    diccionario_global_FFL_I1[f"Coeficiente_Hill_{Hill_modulo_Y}"] = [matrix_info_YX_I1_GLOBAL,matrix_info_ZX_I1_GLOBAL,matrix_info_ZY_I1_GLOBAL,matrix_info_ZXY_I1_GLOBAL]
    np.save('MATRICES_INFORMACION_FFL_I1_AND_Modificacion_Hill_Y_inicial.npy', diccionario_global_FFL_I1)
    diccionario_global_FFL_I3[f"Coeficiente_Hill_{Hill_modulo_Y}"] = [matrix_info_YX_I3_GLOBAL,matrix_info_ZX_I3_GLOBAL,matrix_info_ZY_I3_GLOBAL,matrix_info_ZXY_I3_GLOBAL]
    np.save('MATRICES_INFORMACION_FFL_I3_AND_Modificacion_Hill_Y_inicial.npy', diccionario_global_FFL_I3)

# %%
