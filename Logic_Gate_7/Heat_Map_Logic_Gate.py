#%%
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
from tqdm import tqdm
@njit

def logic_gate_function_7(A, B, K, a):
    term1 = (K[0] * A) / ((K[3] + A) * (K[4] + B))
    term2 = (K[1] * B) / ((K[3] + A) * (K[4] + B))
    term3 = (K[2] * A * B) / ((K[3] + A) * (K[4] + B))

    result = a*(term1 + term2 - term3)

    return result

# Definir valores de K arbitrarios para la simulación
K = [1, 1, 1, 1, 1]

# Crear una malla de valores de A y R
A_values = np.linspace(0, 1, 1000)
R_values = np.linspace(0, 1, 1000)

# Inicializar la matriz de resultados
results = np.zeros((len(A_values), len(R_values)))

# Rellenar la matriz con los resultados de la función
for i, A in tqdm(enumerate(A_values)):
    for j, R in enumerate(R_values):
        results[i, j] = logic_gate_function_7(A, R, K, a = 1)
#%%
# Generar el mapa de calor
plt.figure(figsize=(8, 5))
plt.imshow(results, extent=[0, 1, 0, 1], origin='lower', aspect='auto', cmap='turbo')
plt.colorbar(label='Function Output')
plt.title('Heatmap of Logic Gate 7', fontsize = 15)
plt.xlabel('Y Value', fontsize = 15)
plt.ylabel('X Value', fontsize = 15)
plt.savefig("HeatMap_Logic_Gate_7.png", dpi = 1000)
# %%
