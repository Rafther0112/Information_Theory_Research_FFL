#%%
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
from tqdm import tqdm

@njit
def function_5(X, Y, kxz, Kxz, kyz, Kyz, n):
    term1 = kxz * (X**n) / (X**n + Kxz**n)
    term2 = kyz * (Y**n) / (Y**n + Kyz**n)
    result = term1 + term2
    return 1 - result

# Parámetros arbitrarios
kxz = 0.5
kyz = 0.5
Kxz = 0.5
Kyz = 0.5
n = 2

# Crear una malla de valores de X y Y
X_values = np.linspace(0, 1, 1000)
Y_values = np.linspace(0, 1, 1000)

# Inicializar la matriz de resultados
results = np.zeros((len(X_values), len(Y_values)))

# Rellenar la matriz con los resultados de la función
for i, X in tqdm(enumerate(X_values)):
    for j, Y in enumerate(Y_values):
        results[i, j] = function_5(X, Y, kxz, Kxz, kyz, Kyz, n)
#%%
# Generar el mapa de calor
plt.figure(figsize=(8, 5))
plt.imshow(results, extent=[0, 1, 0, 1], origin='lower', aspect='auto', cmap='turbo')
plt.colorbar(label='Function Output')
plt.title('Heatmap of Logic Gate 5', fontsize = 15)
plt.xlabel('Y Value', fontsize = 15)
plt.ylabel('X Value', fontsize = 15)
plt.savefig("HeatMap_Logic_Gate_5.png", dpi = 1000)


# %%