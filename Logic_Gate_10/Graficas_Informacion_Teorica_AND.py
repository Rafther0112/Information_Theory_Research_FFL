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
#%% Cargamos los datos teoricos encontrados con mathematica para los coherentes
datos_teoricos_FFL_C1 = np.empty(shape=(4, 91), dtype='object')
datos_teoricos_FFL_C2 = np.empty(shape=(4, 91), dtype='object')
datos_teoricos_FFL_C3 = np.empty(shape=(4, 91), dtype='object')
datos_teoricos_FFL_C4 = np.empty(shape=(4, 91), dtype='object')
for position_Hill, hill in enumerate([1,2,3,4]):
    with open(f'Theorical_Results_Information/datos_FFL_C1_Hill{hill}.csv', 'r') as archivo:
        lector = csv.reader(archivo)
        datos_FFL_C1 = list(lector)
        for posicion_Valor, valor in enumerate(datos_FFL_C1[0]):
            dato = float(valor)
            datos_teoricos_FFL_C1[position_Hill][posicion_Valor] = dato

    with open(f'Theorical_Results_Information/datos_FFL_C2_Hill{hill}.csv', 'r') as archivo:
        lector = csv.reader(archivo)
        datos_FFL_C2 = list(lector)
        for posicion_Valor, valor in enumerate(datos_FFL_C2[0]):
            dato = float(valor)
            datos_teoricos_FFL_C2[position_Hill][posicion_Valor] = dato

    with open(f'Theorical_Results_Information/datos_FFL_C3_Hill{hill}.csv', 'r') as archivo:
        lector = csv.reader(archivo)
        datos_FFL_C3 = list(lector)
        for posicion_Valor, valor in enumerate(datos_FFL_C3[0]):
            dato = float(valor)
            datos_teoricos_FFL_C3[position_Hill][posicion_Valor] = dato

    with open(f'Theorical_Results_Information/datos_FFL_C4_Hill{hill}.csv', 'r') as archivo:
        lector = csv.reader(archivo)
        datos_FFL_C4 = list(lector)
        for posicion_Valor, valor in enumerate(datos_FFL_C4[0]):
            dato = float(valor)
            datos_teoricos_FFL_C4[position_Hill][posicion_Valor] = dato
#%% Cargamos los datos teoricos encontrados con mathematica para los incoherentes
datos_teoricos_FFL_I1 = np.empty(shape=(4, 91), dtype='object')
datos_teoricos_FFL_I2 = np.empty(shape=(4, 91), dtype='object')
datos_teoricos_FFL_I3 = np.empty(shape=(4, 91), dtype='object')
datos_teoricos_FFL_I4 = np.empty(shape=(4, 91), dtype='object')
for position_Hill, hill in enumerate([1,2,3,4]):
    with open(f'Theorical_Results_Information/datos_FFL_I1_Hill{hill}.csv', 'r') as archivo:
        lector = csv.reader(archivo)
        datos_FFL_I1 = list(lector)
        for posicion_Valor, valor in enumerate(datos_FFL_I1[0]):
            dato = float(valor)
            datos_teoricos_FFL_I1[position_Hill][posicion_Valor] = dato

    with open(f'Theorical_Results_Information/datos_FFL_I2_Hill{hill}.csv', 'r') as archivo:
        lector = csv.reader(archivo)
        datos_FFL_I2 = list(lector)
        for posicion_Valor, valor in enumerate(datos_FFL_I2[0]):
            dato = float(valor)
            datos_teoricos_FFL_I2[position_Hill][posicion_Valor] = dato

    with open(f'Theorical_Results_Information/datos_FFL_I3_Hill{hill}.csv', 'r') as archivo:
        lector = csv.reader(archivo)
        datos_FFL_I3 = list(lector)
        for posicion_Valor, valor in enumerate(datos_FFL_I3[0]):
            dato = float(valor)
            datos_teoricos_FFL_I3[position_Hill][posicion_Valor] = dato

    with open(f'Theorical_Results_Information/datos_FFL_I4_Hill{hill}.csv', 'r') as archivo:
        lector = csv.reader(archivo)
        datos_FFL_I4 = list(lector)
        for posicion_Valor, valor in enumerate(datos_FFL_I4[0]):
            dato = float(valor)
            datos_teoricos_FFL_I4[position_Hill][posicion_Valor] = dato
# %%
# _________________ANALISIS DE INFORMACION EN ESTADO ESTACIONARIO____________________________
#%% Analisis para los coherentes
sampling = 150
informacion_estado_estacionario_C_X_Z = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]

for posicion_Hill, Valor_Hill in enumerate(range(1,5)):
    for valor_K in range(0,10):

        data_C1 = {'X': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_C1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C1 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C1)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][0].append((1/2)*np.log2((Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2])/(Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2] - (Cov_matrix_C1[0][2])**2)))

        data_C2 = {'X': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_C2[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C2 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C2)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][1].append((1/2)*np.log2((Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2])/(Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2] - (Cov_matrix_C2[0][2])**2)))

        data_C3 = {'X': simulacion_FFL_C3[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_C3[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_C3[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C3 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C3)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][2].append((1/2)*np.log2((Cov_matrix_C3[0][0]* Cov_matrix_C3[2][2])/(Cov_matrix_C3[0][0]* Cov_matrix_C3[2][2] - (Cov_matrix_C3[0][2])**2)))

        data_C4 = {'X': simulacion_FFL_C4[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_C4[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_C4[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C4 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C4)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][3].append((1/2)*np.log2((Cov_matrix_C4[0][0]* Cov_matrix_C4[2][2])/(Cov_matrix_C4[0][0]* Cov_matrix_C4[2][2] - (Cov_matrix_C4[0][2])**2)))


# %% Graficas multiples de informacion coherentes I(X;Z)
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
fig.suptitle(r"Información mutua en estado estacionario $I(\langle X \rangle;\langle Z \rangle)$" + "\n" +  r"  FFL Coherentes con diferente coeficiente de Hill |  Compuerta lógica AND", fontsize = 16)

regulacion_K = 0
# Puedes acceder a cada eje individualmente utilizando la notación de índices
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "C1", color = Red_color,alpha = 0.7)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "C2", color = Orange_color,alpha = 0.7)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "C3", color = Green_color,alpha = 0.7)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "C4", color = Blue_color,alpha = 0.7)

axes[0, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C1[regulacion_K]), color = Red_color, linestyle = "-.")
axes[0, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C2[regulacion_K]), color = Orange_color, linestyle = "-.")
axes[0, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C3[regulacion_K]), color = Green_color, linestyle = "--")
axes[0, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C4[regulacion_K]), color = Blue_color, linestyle = "--")

axes[0, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[0, 0].set_ylim([-0.001, 0.065])
axes[0, 0].legend() 

regulacion_K = 1
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "C1", color = Red_color,alpha = 0.7)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "C2", color = Orange_color,alpha = 0.7)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "C3", color = Green_color,alpha = 0.7)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "C4", color = Blue_color,alpha = 0.7)

axes[0, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C1[regulacion_K]), color = Red_color, linestyle = "-.")
axes[0, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C2[regulacion_K]), color = Orange_color, linestyle = "-.")
axes[0, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C3[regulacion_K]), color = Green_color, linestyle = "--")
axes[0, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C4[regulacion_K]), color = Blue_color, linestyle = "--")

axes[0, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[0, 1].set_ylim([-0.001, 0.12])
axes[0, 1].legend()

regulacion_K = 2
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "C1", color = Red_color,alpha = 0.7)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "C2", color = Orange_color,alpha = 0.7)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "C3", color = Green_color,alpha = 0.7)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "C4", color = Blue_color,alpha = 0.7)

axes[1, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C1[regulacion_K]), color = Red_color, linestyle = "-.")
axes[1, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C2[regulacion_K]), color = Orange_color, linestyle = "-.")
axes[1, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C3[regulacion_K]), color = Green_color, linestyle = "--")
axes[1, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C4[regulacion_K]), color = Blue_color, linestyle = "--")

axes[1, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[1, 0].set_ylim([-0.001, 0.15])
axes[1, 0].legend()

regulacion_K = 3
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "C1", color = Red_color,alpha = 0.7)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "C2", color = Orange_color,alpha = 0.7)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "C3", color = Green_color,alpha = 0.7)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "C4", color = Blue_color,alpha = 0.7)

axes[1, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C1[regulacion_K]), color = Red_color, linestyle = "-.")
axes[1, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C2[regulacion_K]), color = Orange_color, linestyle = "-.")
axes[1, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C3[regulacion_K]), color = Green_color, linestyle = "--")
axes[1, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_C4[regulacion_K]), color = Blue_color, linestyle = "--")

axes[1, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[1, 1].set_ylim([-0.001, 0.165]) 
axes[1, 1].legend()

plt.tight_layout()
plt.savefig("Imagenes_Resultados/Informacion_Estacionaria_Coherente_Multiple_I(X_Z)_AND.jpg", dpi=1000)

#%% Analisis para los incoherentes
sampling = 150
informacion_estado_estacionario_C_X_Z = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]
informacion_estado_estacionario_C_Z_X_Y = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]

for posicion_Hill, Valor_Hill in enumerate(range(1,5)):
    for valor_K in range(0,10):

        data_C1 = {'X': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_I1[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C1 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C1)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][0].append((1/2)*np.log2((Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2])/(Cov_matrix_C1[0][0]* Cov_matrix_C1[2][2] - (Cov_matrix_C1[0][2])**2)))

        data_C2 = {'X': simulacion_FFL_I2[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_I2[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_I2[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C2 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C2)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][1].append((1/2)*np.log2((Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2])/(Cov_matrix_C2[0][0]* Cov_matrix_C2[2][2] - (Cov_matrix_C2[0][2])**2)))

        data_C3 = {'X': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_I3[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C3 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C3)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][2].append((1/2)*np.log2((Cov_matrix_C3[0][0]* Cov_matrix_C3[2][2])/(Cov_matrix_C3[0][0]* Cov_matrix_C3[2][2] - (Cov_matrix_C3[0][2])**2)))

        data_C4 = {'X': simulacion_FFL_I4[f"Coeficiente_Hill_{Valor_Hill}"][0][valor_K][:,sampling:].flatten(),
                   'Y': simulacion_FFL_I4[f"Coeficiente_Hill_{Valor_Hill}"][1][valor_K][:,sampling:].flatten(),
                   'Z': simulacion_FFL_I4[f"Coeficiente_Hill_{Valor_Hill}"][2][valor_K][:,sampling:].flatten()}
        Cov_matrix_C4 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C4)))
        informacion_estado_estacionario_C_X_Z[posicion_Hill][3].append((1/2)*np.log2((Cov_matrix_C4[0][0]* Cov_matrix_C4[2][2])/(Cov_matrix_C4[0][0]* Cov_matrix_C4[2][2] - (Cov_matrix_C4[0][2])**2)))
# %% Graficas multiples de informacion Incoherentes I(X;Z)

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
fig.suptitle(r"Información mutua en estado estacionario $I(\langle X \rangle;\langle Z \rangle)$" + "\n" +  r"  FFL Incoherentes con diferente coeficiente de Hill |  Compuerta lógica AND", fontsize = 16)

regulacion_K = 0
# Puedes acceder a cada eje individualmente utilizando la notación de índices
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "I1", color = Red_color,alpha = 0.7)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "I2", color = Orange_color,alpha = 0.7)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "I3", color = Green_color,alpha = 0.7)
axes[0, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "I4", color = Blue_color,alpha = 0.7)

axes[0, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I1[regulacion_K]), color = Red_color, linestyle = "--")
axes[0, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I2[regulacion_K]), color = Orange_color, linestyle = "-.")
axes[0, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I3[regulacion_K]), color = Green_color, linestyle = "--")
axes[0, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I4[regulacion_K]), color = Blue_color, linestyle = "-.")

axes[0, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[0, 0].set_ylim([-0.005, 0.06]) 
axes[0, 0].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[0, 0].legend() 

regulacion_K = 1
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "I1", color = Red_color,alpha = 0.7)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "I2", color = Orange_color,alpha = 0.7)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "I3", color = Green_color,alpha = 0.7)
axes[0, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "I4", color = Blue_color,alpha = 0.7)

axes[0, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I1[regulacion_K]), color = Red_color, linestyle = "--")
axes[0, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I2[regulacion_K]), color = Orange_color, linestyle = "-.")
axes[0, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I3[regulacion_K]), color = Green_color, linestyle = "--")
axes[0, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I4[regulacion_K]), color = Blue_color, linestyle = "-.")

axes[0, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[0, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[0, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[0, 1].set_ylim([-0.005, 0.135]) 
axes[0, 1].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[0, 1].legend()

regulacion_K = 2
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "I1", color = Red_color,alpha = 0.7)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "I2", color = Orange_color,alpha = 0.7)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "I3", color = Green_color,alpha = 0.7)
axes[1, 0].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "I4", color = Blue_color,alpha = 0.7)

axes[1, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I1[regulacion_K]), color = Red_color, linestyle = "--")
axes[1, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I2[regulacion_K]), color = Orange_color, linestyle = "-.")
axes[1, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I3[regulacion_K]), color = Green_color, linestyle = "--")
axes[1, 0].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I4[regulacion_K]), color = Blue_color, linestyle = "-.")

axes[1, 0].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 0].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 0].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[1, 0].set_ylim([-0.005, 0.2]) 
axes[1, 0].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[1, 0].legend()

regulacion_K = 3
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][0]), label = "I1", color = Red_color,alpha = 0.7)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][1]), label = "I2", color = Orange_color,alpha = 0.7)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][2]), label = "I3", color = Green_color,alpha = 0.7)
axes[1, 1].scatter(range(1,11), np.array(informacion_estado_estacionario_C_X_Z[regulacion_K][3]), label = "I4", color = Blue_color,alpha = 0.7)

axes[1, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I1[regulacion_K]), color = Red_color, linestyle = "--")
axes[1, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I2[regulacion_K]), color = Orange_color, linestyle = "-.")
axes[1, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I3[regulacion_K]), color = Green_color, linestyle = "--")
axes[1, 1].plot(np.arange(1,10.1,0.1), np.array(datos_teoricos_FFL_I4[regulacion_K]), color = Blue_color, linestyle = "-.")

axes[1, 1].set_title(rf'Coeficiente de Hill n = {regulacion_K+1}', fontsize = 14)
axes[1, 1].set_xlabel(r"Valor tasa $K_x$", fontsize = 14)
axes[1, 1].set_ylabel(r"$I(\langle X \rangle; \langle Z \rangle)$", fontsize = 14)
axes[1, 1].set_ylim([-0.005, 0.22]) 
axes[1, 1].axhline(y = 0, linestyle = "--", color = "gray", label = "Limite \n inferior")
axes[1, 1].legend()

plt.tight_layout()
plt.savefig("Imagenes_Resultados/Informacion_Estacionaria_Incoherente_Multiple_I(X_Z)_AND.jpg", dpi=1000)
#%%



