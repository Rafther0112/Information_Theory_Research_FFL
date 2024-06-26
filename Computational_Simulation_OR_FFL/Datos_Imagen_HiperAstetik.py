#%% IMPORTS Y FUNCIONES
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

# %% ANALISIS TEMPORAL PARA COHERENTES I(X:Z|Y)
tiempo_propio_X = np.arange(0,350, 2)
tiempo_propio_Y = np.arange(0,350, 2)
tiempo_propio_Z = np.arange(0,350, 2)


informacion_configuracion_X_ZY_C1 = np.empty(shape=(len(tiempo_propio_X),len(tiempo_propio_Y),len(tiempo_propio_Z)), dtype='float')


for posicionX, TauX in enumerate((tqdm(tiempo_propio_X))):
    for posicionY, TauY in enumerate(((tiempo_propio_Y))):
        for posicionZ, TauZ in enumerate(tiempo_propio_Z):
            data_C1 = {'X': simulacion_FFL_C1[f"Coeficiente_Hill_{4}"][0][3][:,TauX],
                    'Y': simulacion_FFL_C1[f"Coeficiente_Hill_{4}"][1][3][:,TauY],
                    'Z': simulacion_FFL_C1[f"Coeficiente_Hill_{4}"][2][3][:,TauZ]}
            Cov_matrix_C1 = np.array(pd.DataFrame.cov(pd.DataFrame(data_C1)))
            informacion_configuracion_X_ZY_C1[posicionX][posicionY][posicionZ] = (1/2)*np.log2((Cov_matrix_C1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C1))/((Cov_matrix_C1[0][0]*conditional_covarianza_Z_Y(Cov_matrix_C1))-(conditional_covarianza_X_ZdadoY(Cov_matrix_C1))**2))


#%%

# Crear un array 3D con valores según una función Gaussiana centrada en la diagonal
array_3d = informacion_configuracion_X_ZY_C1

# Obtener las coordenadas para los puntos del array
x, y, z = np.meshgrid(range(array_3d.shape[0]), range(array_3d.shape[1]), range(array_3d.shape[2]))

# Crear la figura 3D
fig = plt.figure(figsize=(10,6))

# Crear un objeto de la clase Axes3D
ax = fig.add_subplot(111, projection='3d')

# Graficar los puntos con colores según los valores en el array y usar cmap 'viridis'
scatter = ax.scatter(y, x, z, c=array_3d.flatten(), cmap='Spectral')

# Ajustar las etiquetas de los ejes
ax.set_xlabel(r" Tiempo propio $\tau_Y$")
ax.set_ylabel(r" Tiempo propio $\tau_X$")
ax.set_zlabel(r" Tiempo propio $\tau_Z$")

ax.set_title(r"Información condicional temporal $I(X;Z|Y)$")
# Agregar la barra de colores
cbar = plt.colorbar(scatter, label= r"Valor de información $I(X;Z|Y)$")

plt.tight_layout()
plt.savefig("Informacion_Temporal_I_XZ|Y_3D.jpg", dpi = 1000)
#plt.savefig("Informacion_Temporal.jpg", dpi = 1000)


# %%
plt.imshow(array_3d[40])
# %%

# %%
