#%% IMPORTS Y FUNCIONES
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
from tqdm import tqdm
Red_color = "#d62728"
Orange_color = "#ff7f0e"
Green_color = "#2ca02c"
Blue_color =  "#1f77b4"

def conditional_covarianza_Z_Y(covariance_matrix):
    return covariance_matrix[2][2] - ((covariance_matrix[2][1])**2)/(covariance_matrix[1][1])

def conditional_covarianza_X_ZdadoY(covariance_matrix):
    return covariance_matrix[0][2] - ((covariance_matrix[0][1])*(covariance_matrix[2][1]))/(covariance_matrix[1][1])

# %% IMPORTAMOS DATOS COMPUTACIONALES
simulacion_FFL_C1 = np.load('Resultados_Simulacion/Simulacion_FFL_C1_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_C2 = np.load('Resultados_Simulacion/Simulacion_FFL_C2_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_C3 = np.load('Resultados_Simulacion/Simulacion_FFL_C3_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_C4 = np.load('Resultados_Simulacion/Simulacion_FFL_C4_OR_final.npy', allow_pickle=True).item()

simulacion_FFL_I1 = np.load('Resultados_Simulacion/Simulacion_FFL_I1_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_I2 = np.load('Resultados_Simulacion/Simulacion_FFL_I2_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_I3 = np.load('Resultados_Simulacion/Simulacion_FFL_I3_OR_final.npy', allow_pickle=True).item()
simulacion_FFL_I4 = np.load('Resultados_Simulacion/Simulacion_FFL_I4_OR_final.npy', allow_pickle=True).item()
# %% ANALISIS TEMPORAL PARA COHERENTES I(X;Z)
tiempo_propio_X = np.arange(0,350, 2)
tiempo_propio_Y = np.arange(0,350, 2)
tiempo_propio_Z = np.arange(0,350, 2)

informacion_promedio_Maxima_Hill_C1 = np.empty(shape=(4,10), dtype='float')
informacion_promedio_Maxima_Hill_C2 = np.empty(shape=(4,10), dtype='float')

for posicion_Hill, Valor_Hill in enumerate((range(1,5))):
    informacion_promedio_Maxima_C1 = []
    informacion_promedio_Maxima_C2 = []
    for posicion_K, valor_K in enumerate((tqdm(range(0,10)))):
        #Cuando entramos aqui estamos fijando una combinacion (Kx,Hill) es decir, un estado especifico de la
        #configuracion. Lo que hacemos hacia abajo solamente es la optimizacion temporal
        informacion_configuracion_X_Z_C1 = np.empty(shape=(len(tiempo_propio_Z),len(tiempo_propio_X)), dtype='float')
        informacion_configuracion_X_ZY_C1 = np.empty(shape=(len(tiempo_propio_Z),len(tiempo_propio_X)), dtype='float')

        informacion_configuracion_X_Z_C2 = np.empty(shape=(len(tiempo_propio_Z),len(tiempo_propio_X)), dtype='float')
        informacion_configuracion_X_ZY_C2 = np.empty(shape=(len(tiempo_propio_Z),len(tiempo_propio_X)), dtype='float')
    
        for posicionX, TauX in enumerate(((tiempo_propio_X))):
            for posicionZ, TauZ in enumerate(tiempo_propio_Z):
                data_C1 = {'X': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,TauX],
                        'Y': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,-1],
                        'Z': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,TauZ]}
                Cov_matrix_C1 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C1)))
                informacion_configuracion_X_Z_C1[posicionX][posicionZ] = (1/2)*np.log2((Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2])/(Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2] - (Cov_matrix_C1[0][2])**2))
                informacion_configuracion_X_ZY_C1[posicionX][posicionZ] = (1/2)*np.log2((Cov_matrix_C1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C1))/((Cov_matrix_C1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C1))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C1))**2))

                data_C2 = {'X': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,TauX],
                        'Y': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,-1],
                        'Z': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,TauZ]}
                Cov_matrix_C2 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C2)))
                informacion_configuracion_X_Z_C2[posicionX][posicionZ] = (1/2)*np.log2((Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2])/(Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2] - (Cov_matrix_C2[0][2])**2))
                informacion_configuracion_X_ZY_C2[posicionX][posicionZ] = (1/2)*np.log2((Cov_matrix_C2[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C2))/((Cov_matrix_C2[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C2))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C2))**2))
        
        informacion_maxima_XZ_C1 = []
        tiempos_maximos_XZ_C1   = []
        data_C1 = informacion_configuracion_X_Z_C1
        # Itera sobre cada columna (valor de X)
        for x in range(data_C1.shape[1]):
            # Obtén la columna actual
            columna_x = data_C1[:, x]

            if np.all(np.isnan(columna_x)):
                valor_maximo = 0
                tiempo_maximo = (x, x)
            else:
                indice_maximo = np.nanargmax(columna_x)
                valor_maximo = columna_x[indice_maximo]
                tiempo_maximo = (x, indice_maximo)

            informacion_maxima_XZ_C1.append(valor_maximo)
            tiempos_maximos_XZ_C1.append(tiempo_maximo)
        informacion_promedio_Maxima_Hill_C1[posicion_Hill][posicion_K] = np.mean(informacion_maxima_XZ_C1)

        informacion_maxima_XZ_C2 = []
        tiempos_maximos_XZ_C2   = []
        data_C2 = informacion_configuracion_X_Z_C2
        # Itera sobre cada columna (valor de X)
        for x in range(data_C2.shape[1]):
            # Obtén la columna actual
            columna_x = data_C2[:, x]

            if np.all(np.isnan(columna_x)):
                valor_maximo = 0
                tiempo_maximo = (x, x)
            else:
                indice_maximo = np.nanargmax(columna_x)
                valor_maximo = columna_x[indice_maximo]
                tiempo_maximo = (x, indice_maximo)

            informacion_maxima_XZ_C2.append(valor_maximo)
            tiempos_maximos_XZ_C2.append(tiempo_maximo)
        informacion_promedio_Maxima_Hill_C2[posicion_Hill][posicion_K] = np.mean(informacion_maxima_XZ_C2)
#%% GRAFICAS TEMPORAL COHERENTES I(X;Z)
plt.figure(figsize=(8,3))
plt.xlabel(r'Valor tasa $K_X$')
plt.ylabel(r'Coeficiente de Hill n')
plt.title(r'Información maxima temporal ($C_1$ , $C_4$) $I_{max}(X;Z)$' + "\n" + r"Compuerta lógica OR")

indices_filas = np.arange(1, informacion_promedio_Maxima_Hill_C1.shape[0] + 1)
indices_columnas = np.arange(1, informacion_promedio_Maxima_Hill_C1.shape[1] + 1)

plt.imshow(informacion_promedio_Maxima_Hill_C1, cmap='Spectral')
plt.xticks(np.arange(len(indices_columnas)), indices_columnas)
plt.yticks(np.arange(len(indices_filas)), indices_filas)
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.set_label(r'Valor de Información')
#plt.savefig("Imagenes_Resultados/Informacion_maxima_temporal_C1C4_color.jpg", dpi = 1000)
np.save('Resultados_Matrices_Color/Informacion_maxima_temporal_IXZ_C1C4_OR_color.npy', informacion_promedio_Maxima_Hill_C1)
#%%
informacion_promedio_Maxima_Hill_C1
#%%

plt.figure(figsize=(8,3))
plt.xlabel(r'Valor tasa $K_X$')
plt.ylabel(r'Coeficiente de Hill n')
plt.title(r'Información maxima temporal ($C_2$ , $C_3$) $I_{max}(X;Z)$'  + "\n" + r"Compuerta lógica OR")

indices_filas = np.arange(1, informacion_promedio_Maxima_Hill_C2.shape[0] + 1)
indices_columnas = np.arange(1, informacion_promedio_Maxima_Hill_C2.shape[1] + 1)

plt.imshow(informacion_promedio_Maxima_Hill_C2, cmap='Spectral')
plt.xticks(np.arange(len(indices_columnas)), indices_columnas)
plt.yticks(np.arange(len(indices_filas)), indices_filas)
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.set_label(r'Valor de Información')
#plt.savefig("Imagenes_Resultados/Informacion_maxima_temporal_C2C3_color.jpg", dpi = 1000)
np.save('Resultados_Matrices_Color/Informacion_maxima_temporal_IXZ_C2C3_OR_color.npy', informacion_promedio_Maxima_Hill_C2)
#%%
informacion_promedio_Maxima_Hill_C2
#%%

#%% ANALISIS TEMPORAL PARA INCOHERENTES I(X;Z)

tiempo_propio_X = np.arange(0,350, 2)
tiempo_propio_Y = np.arange(0,350, 2)
tiempo_propio_Z = np.arange(0,350, 2)

informacion_promedio_Maxima_Hill_I1 = np.empty(shape=(4,10), dtype='float')
informacion_promedio_Maxima_Hill_I3 = np.empty(shape=(4,10), dtype='float')

for posicion_Hill, Valor_Hill in enumerate((range(1,5))):
    informacion_promedio_Maxima_I1 = []
    informacion_promedio_Maxima_I3 = []
    for posicion_K, valor_K in enumerate((tqdm(range(0,10)))):
        #Cuando entramos aqui estamos fijando una combinacion (Kx,Hill) es decir, un estado especifico de la
        #configuracion. Lo que hacemos hacia abajo solamente es la optimizacion temporal
        informacion_configuracion_X_Z_I1 = np.empty(shape=(len(tiempo_propio_Z),len(tiempo_propio_X)), dtype='float')
        informacion_configuracion_X_ZY_I1 = np.empty(shape=(len(tiempo_propio_Z),len(tiempo_propio_X)), dtype='float')

        informacion_configuracion_X_Z_I3 = np.empty(shape=(len(tiempo_propio_Z),len(tiempo_propio_X)), dtype='float')
        informacion_configuracion_X_ZY_I3 = np.empty(shape=(len(tiempo_propio_Z),len(tiempo_propio_X)), dtype='float')
    
        for posicionX, TauX in enumerate(((tiempo_propio_X))):
            for posicionZ, TauZ in enumerate(tiempo_propio_Z):
                data_I1 = {'X': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,TauX],
                        'Y': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,-1],
                        'Z': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,TauZ]}
                Cov_matrix_I1 = np.array(pd.DataFrame.cov(pd.DataFrame(data_I1)))
                informacion_configuracion_X_Z_I1[posicionX][posicionZ] = (1/2)*np.log2((Cov_matrix_I1[0][0]* Cov_matrix_I1[2][2])/(Cov_matrix_I1[0][0]* Cov_matrix_I1[2][2] - (Cov_matrix_I1[0][2])**2))
                informacion_configuracion_X_ZY_I1[posicionX][posicionZ] = (1/2)*np.log2((Cov_matrix_I1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_I1))/((Cov_matrix_I1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_I1))-(conditional_covarianza_X_ZdadoY(Cov_matrix_I1))**2))

                data_I3 = {'X': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,TauX],
                        'Y': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,-1],
                        'Z': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,TauZ]}
                Cov_matrix_I3 = np.array(pd.DataFrame.cov(pd.DataFrame(data_I3)))
                informacion_configuracion_X_Z_I3[posicionX][posicionZ] = (1/2)*np.log2((Cov_matrix_I3[0][0]* Cov_matrix_I3[2][2])/(Cov_matrix_I3[0][0]* Cov_matrix_I3[2][2] - (Cov_matrix_I3[0][2])**2))
                informacion_configuracion_X_ZY_I3[posicionX][posicionZ] = (1/2)*np.log2((Cov_matrix_I3[0][0]*conditional_covarianza_Z_Y(Cov_matrix_I3))/((Cov_matrix_I3[0][0]*conditional_covarianza_Z_Y(Cov_matrix_I3))-(conditional_covarianza_X_ZdadoY(Cov_matrix_I3))**2))

        informacion_maxima_XZ_I1 = []
        tiempos_maximos_XZ_I1   = []
        data_I1 = informacion_configuracion_X_Z_I1
        # Itera sobre cada columna (valor de X)
        for x in range(data_I1.shape[1]):
            # Obtén la columna actual
            columna_x = data_I1[:, x]

            if np.all(np.isnan(columna_x)):
                valor_maximo = 0
                tiempo_maximo = (x, x)
            else:
                indice_maximo = np.nanargmax(columna_x)
                valor_maximo = columna_x[indice_maximo]
                tiempo_maximo = (x, indice_maximo)

            informacion_maxima_XZ_I1.append(valor_maximo)
            tiempos_maximos_XZ_I1.append(tiempo_maximo)
        informacion_promedio_Maxima_Hill_I1[posicion_Hill][posicion_K] = np.mean(informacion_maxima_XZ_I1)

        informacion_maxima_XZ_I3 = []
        tiempos_maximos_XZ_I3   = []
        data_I3 = informacion_configuracion_X_Z_I3
        # Itera sobre cada columna (valor de X)
        for x in range(data_I3.shape[1]):
            # Obtén la columna actual
            columna_x = data_I3[:, x]

            if np.all(np.isnan(columna_x)):
                valor_maximo = 0
                tiempo_maximo = (x, x)
            else:
                indice_maximo = np.nanargmax(columna_x)
                valor_maximo = columna_x[indice_maximo]
                tiempo_maximo = (x, indice_maximo)

            informacion_maxima_XZ_I3.append(valor_maximo)
            tiempos_maximos_XZ_I3.append(tiempo_maximo)
        informacion_promedio_Maxima_Hill_I3[posicion_Hill][posicion_K] = np.mean(informacion_maxima_XZ_I3)
# %% GRAFICA TEMPORAL INCOHERENTES I(X;Z)
plt.figure(figsize=(8,3))
plt.xlabel(r'Valor tasa $K_X$')
plt.ylabel(r'Coeficiente de Hill n')
plt.title(r'Información maxima temporal ($I_1$ , $I_2$) $I_{max}(X;Z)$'  + "\n" + r"Compuerta lógica OR")
indices_filas = np.arange(1, informacion_promedio_Maxima_Hill_I1.shape[0] + 1)
indices_columnas = np.arange(1, informacion_promedio_Maxima_Hill_I1.shape[1] + 1)
plt.imshow(informacion_promedio_Maxima_Hill_I1, cmap='Spectral')
plt.xticks(np.arange(len(indices_columnas)), indices_columnas)
plt.yticks(np.arange(len(indices_filas)), indices_filas)
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.set_label(r'Valor de Información')
#plt.savefig("Imagenes_Resultados/Informacion_maxima_temporal_I1I2_color.jpg", dpi = 1000)
np.save('Resultados_Matrices_Color/Informacion_maxima_temporal_IXZ_I1I2_OR_color.npy', informacion_promedio_Maxima_Hill_I1)
#%%
informacion_promedio_Maxima_Hill_I1
#%%

plt.figure(figsize=(8,3))
plt.xlabel(r'Valor tasa $K_X$')
plt.ylabel(r'Coeficiente de Hill n')
plt.title(r'Información maxima temporal ($I_3$ , $I_4$) $I_{max}(X;Z)$'  + "\n" + r"Compuerta lógica OR")
indices_filas = np.arange(1, informacion_promedio_Maxima_Hill_I3.shape[0] + 1)
indices_columnas = np.arange(1, informacion_promedio_Maxima_Hill_I3.shape[1] + 1)
plt.imshow(informacion_promedio_Maxima_Hill_I3, cmap='Spectral')
plt.xticks(np.arange(len(indices_columnas)), indices_columnas)
plt.yticks(np.arange(len(indices_filas)), indices_filas)
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.set_label(r'Valor de Información')
#plt.savefig("Imagenes_Resultados/Informacion_maxima_temporal_I3I4_color.jpg", dpi = 1000)
np.save('Resultados_Matrices_Color/Informacion_maxima_temporal_IXZ_I3I4_OR_color.npy', informacion_promedio_Maxima_Hill_I3)
#%%
informacion_promedio_Maxima_Hill_I3
#%%
# %% ANALISIS TEMPORAL PARA COHERENTES I(X:Z|Y)
tiempo_propio_X = np.arange(0,350, 100)
tiempo_propio_Y = np.arange(0,350, 100)
tiempo_propio_Z = np.arange(0,350, 100)

informacion_promedio_Maxima_Hill_C1 = np.empty(shape=(4,10), dtype='float')
informacion_promedio_Maxima_Hill_C2 = np.empty(shape=(4,10), dtype='float')

for posicion_Hill, Valor_Hill in enumerate((range(1,5))):
    informacion_promedio_Maxima_I_X_ZY_C1 = []
    informacion_promedio_Maxima_I_X_ZY_C2 = []
    for posicion_K, valor_K in enumerate((tqdm(range(0,10)))):
        #Cuando entramos aqui estamos fijando una combinacion (Kx,Hill) es decir, un estado especifico de la
        #configuracion. Lo que hacemos hacia abajo solamente es la optimizacion temporal
        informacion_configuracion_X_Z_C1 = np.empty(shape=(len(tiempo_propio_X),len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='float')
        informacion_configuracion_X_ZY_C1 = np.empty(shape=(len(tiempo_propio_X),len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='float')

        informacion_configuracion_X_Z_C2 = np.empty(shape=(len(tiempo_propio_X),len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='float')
        informacion_configuracion_X_ZY_C2 = np.empty(shape=(len(tiempo_propio_X),len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='float')

        for posicionX, TauX in enumerate(((tiempo_propio_X))):
            for posicionY, TauY in enumerate(((tiempo_propio_Y))):
                for posicionZ, TauZ in enumerate(tiempo_propio_Z):
                    data_C1 = {'X': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,TauX],
                            'Y': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,TauY],
                            'Z': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,TauZ]}
                    Cov_matrix_C1 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C1)))
                    informacion_configuracion_X_Z_C1[posicionX][posicionY][posicionZ] = (1/2)*np.log2((Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2])/(Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2] - (Cov_matrix_C1[0][2])**2))
                    informacion_configuracion_X_ZY_C1[posicionX][posicionY][posicionZ] = (1/2)*np.log2((Cov_matrix_C1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C1))/((Cov_matrix_C1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C1))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C1))**2))

                    data_C2 = {'X': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,TauX],
                            'Y': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,-1],
                            'Z': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,TauZ]}
                    Cov_matrix_C2 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C2)))
                    informacion_configuracion_X_Z_C2[posicionX][posicionY][posicionZ] = (1/2)*np.log2((Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2])/(Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2] - (Cov_matrix_C2[0][2])**2))
                    informacion_configuracion_X_ZY_C2[posicionX][posicionY][posicionZ] = (1/2)*np.log2((Cov_matrix_C2[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C2))/((Cov_matrix_C2[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C2))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C2))**2))
        
        maximos_en_z = np.nanmax(informacion_configuracion_X_ZY_C1, axis=2)
        promedio_en_y = np.nanmean(maximos_en_z, axis=1)
        promedio_general_X_ZY_C1 = np.nanmean(promedio_en_y)

        informacion_promedio_Maxima_Hill_C1[posicion_Hill][posicion_K] = promedio_general_X_ZY_C1


        maximos_en_z = np.nanmax(informacion_configuracion_X_ZY_C2, axis=2)
        promedio_en_y = np.nanmean(maximos_en_z, axis=1)
        promedio_general_X_ZY_C2 = np.nanmean(promedio_en_y)

        informacion_promedio_Maxima_Hill_C2[posicion_Hill][posicion_K] = promedio_general_X_ZY_C2
#%%
#%% GRAFICAS TEMPORAL COHERENTES I(X;Z|Y)
plt.figure(figsize=(8,3))
plt.xlabel(r'Valor tasa $K_X$')
plt.ylabel(r'Coeficiente de Hill n')
plt.title(r'Información maxima temporal ($C_1$ , $C_4$) $I_{max}(X;Z|Y)$' + "\n" + r"Compuerta lógica OR")

indices_filas = np.arange(1, informacion_promedio_Maxima_Hill_C1.shape[0] + 1)
indices_columnas = np.arange(1, informacion_promedio_Maxima_Hill_C1.shape[1] + 1)

plt.imshow(informacion_promedio_Maxima_Hill_C1, cmap='Spectral')
plt.xticks(np.arange(len(indices_columnas)), indices_columnas)
plt.yticks(np.arange(len(indices_filas)), indices_filas)
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.set_label(r'Valor de Información')
plt.savefig("Imagenes_Resultados/Informacion_maxima_temporal_I_XZ|Y_C1C4_OR_color.jpg", dpi = 1000)

plt.figure(figsize=(8,3))
plt.xlabel(r'Valor tasa $K_X$')
plt.ylabel(r'Coeficiente de Hill n')
plt.title(r'Información maxima temporal ($C_2$ , $C_3$) $I_{max}(X;Z|Y)$'  + "\n" + r"Compuerta lógica OR")

indices_filas = np.arange(1, informacion_promedio_Maxima_Hill_C2.shape[0] + 1)
indices_columnas = np.arange(1, informacion_promedio_Maxima_Hill_C2.shape[1] + 1)

plt.imshow(informacion_promedio_Maxima_Hill_C2, cmap='Spectral')
plt.xticks(np.arange(len(indices_columnas)), indices_columnas)
plt.yticks(np.arange(len(indices_filas)), indices_filas)
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.set_label(r'Valor de Información')
plt.savefig("Imagenes_Resultados/Informacion_maxima_temporal_I_XZ|Y_C2C3_OR_color.jpg", dpi = 1000)
#%% ANALISIS TEMPORAL PARA INCOHERENTES I(X;Z|Y)
tiempo_propio_X = np.arange(0,350, 25)
tiempo_propio_Y = np.arange(0,350, 25)
tiempo_propio_Z = np.arange(0,350, 25)

informacion_promedio_Maxima_Hill_I1 = np.empty(shape=(4,10), dtype='float')
informacion_promedio_Maxima_Hill_I3 = np.empty(shape=(4,10), dtype='float')

for posicion_Hill, Valor_Hill in enumerate((range(1,5))):
    informacion_promedio_Maxima_I_X_ZY_I1 = []
    informacion_promedio_Maxima_I_X_ZY_I3 = []
    for posicion_K, valor_K in enumerate((tqdm(range(0,10)))):
        #Cuando entramos aqui estamos fijando una combinacion (Kx,Hill) es decir, un estado especifico de la
        #configuracion. Lo que hacemos hacia abajo solamente es la optimizacion temporal
        informacion_configuracion_X_Z_I1 = np.empty(shape=(len(tiempo_propio_X),len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='float')
        informacion_configuracion_X_ZY_I1 = np.empty(shape=(len(tiempo_propio_X),len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='float')

        informacion_configuracion_X_Z_I3 = np.empty(shape=(len(tiempo_propio_X),len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='float')
        informacion_configuracion_X_ZY_I3 = np.empty(shape=(len(tiempo_propio_X),len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='float')

        for posicionX, TauX in enumerate(((tiempo_propio_X))):
            for posicionY, TauY in enumerate(((tiempo_propio_Y))):
                for posicionZ, TauZ in enumerate(tiempo_propio_Z):
                    data_I1 = {'X': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,TauX],
                            'Y': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,TauY],
                            'Z': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,TauZ]}
                    Cov_matrix_I1 = np.array(pd.DataFrame.cov(pd.DataFrame(data_I1)))
                    informacion_configuracion_X_Z_I1[posicionX][posicionY][posicionZ] = (1/2)*np.log2((Cov_matrix_I1[0][0]* Cov_matrix_I1[2][2])/(Cov_matrix_I1[0][0]* Cov_matrix_I1[2][2] - (Cov_matrix_I1[0][2])**2))
                    informacion_configuracion_X_ZY_I1[posicionX][posicionY][posicionZ] = (1/2)*np.log2((Cov_matrix_I1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_I1))/((Cov_matrix_I1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_I1))-(conditional_covarianza_X_ZdadoY(Cov_matrix_I1))**2))

                    data_I3 = {'X': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,TauX],
                            'Y': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,-1],
                            'Z': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,TauZ]}
                    Cov_matrix_I3 = np.array(pd.DataFrame.cov(pd.DataFrame(data_I3)))
                    informacion_configuracion_X_Z_I3[posicionX][posicionY][posicionZ] = (1/2)*np.log2((Cov_matrix_I3[0][0]* Cov_matrix_I3[2][2])/(Cov_matrix_I3[0][0]* Cov_matrix_I3[2][2] - (Cov_matrix_I3[0][2])**2))
                    informacion_configuracion_X_ZY_I3[posicionX][posicionY][posicionZ] = (1/2)*np.log2((Cov_matrix_I3[0][0]*conditional_covarianza_Z_Y(Cov_matrix_I3))/((Cov_matrix_I3[0][0]*conditional_covarianza_Z_Y(Cov_matrix_I3))-(conditional_covarianza_X_ZdadoY(Cov_matrix_I3))**2))
        
        maximos_en_z = np.nanmax(informacion_configuracion_X_ZY_I1, axis=2)
        promedio_en_y = np.nanmean(maximos_en_z, axis=1)
        promedio_general_X_ZY_I1 = np.nanmean(promedio_en_y)

        informacion_promedio_Maxima_Hill_I1[posicion_Hill][posicion_K] = promedio_general_X_ZY_I1


        maximos_en_z = np.nanmax(informacion_configuracion_X_ZY_I3, axis=2)
        promedio_en_y = np.nanmean(maximos_en_z, axis=1)
        promedio_general_X_ZY_I3 = np.nanmean(promedio_en_y)

        informacion_promedio_Maxima_Hill_I3[posicion_Hill][posicion_K] = promedio_general_X_ZY_I3
#%% GRAFICAS TEMPORAL INCOHERENTES I(X;Z|Y)
plt.figure(figsize=(8,3))
plt.xlabel(r'Valor tasa $K_X$')
plt.ylabel(r'Coeficiente de Hill n')
plt.title(r'Información maxima temporal ($I_1$ , $I_2$) $I_{max}(X;Z|Y)$' + "\n" + r"Compuerta lógica OR")

indices_filas = np.arange(1, informacion_promedio_Maxima_Hill_I1.shape[0] + 1)
indices_columnas = np.arange(1, informacion_promedio_Maxima_Hill_I1.shape[1] + 1)

plt.imshow(informacion_promedio_Maxima_Hill_I1, cmap='Spectral')
plt.xticks(np.arange(len(indices_columnas)), indices_columnas)
plt.yticks(np.arange(len(indices_filas)), indices_filas)
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.set_label(r'Valor de Información')
plt.savefig("Imagenes_Resultados/Informacion_maxima_temporal_I_XZ|Y_I1I2_OR_color.jpg", dpi = 1000)

plt.figure(figsize=(8,3))
plt.xlabel(r'Valor tasa $K_X$')
plt.ylabel(r'Coeficiente de Hill n')
plt.title(r'Información maxima temporal ($I_3$ , $I_4$) $I_{max}(X;Z|Y)$'  + "\n" + r"Compuerta lógica OR")

indices_filas = np.arange(1, informacion_promedio_Maxima_Hill_I3.shape[0] + 1)
indices_columnas = np.arange(1, informacion_promedio_Maxima_Hill_I3.shape[1] + 1)

plt.imshow(informacion_promedio_Maxima_Hill_I3, cmap='Spectral')
plt.xticks(np.arange(len(indices_columnas)), indices_columnas)
plt.yticks(np.arange(len(indices_filas)), indices_filas)
plt.gca().invert_yaxis()
cbar = plt.colorbar()
cbar.set_label(r'Valor de Información')
plt.savefig("Imagenes_Resultados/Informacion_maxima_temporal_I_XZ|Y_I3I4_OR_color.jpg", dpi = 1000)

# %%