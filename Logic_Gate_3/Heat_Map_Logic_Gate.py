#%%
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
from tqdm import tqdm
@njit
def logic_gate_function_3(A, R, K):
    # Asegurarse de que K sea un array de numpy para operaciones elemento por elemento
    K = np.array(K)

    term1 = A * K[0] + K[1] + R * K[4]
    term2 = (R * K[2] + K[3]) * (K[1] + R * K[4])
    term3 = A * K[0] * (K[3] + R * K[4])

    numerator1 = term1
    denominator = term2 + term3
    fraction1 = numerator1 / denominator

    numerator2 = K[1] + R * K[4]
    fraction2 = numerator2 / denominator

    result = K[3] * (fraction1 - fraction2)

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
        results[i, j] = logic_gate_function_3(A, R, K)
#%%
# Generar el mapa de calor
plt.figure(figsize=(8, 5))
plt.imshow(results, extent=[0, 1, 0, 1], origin='lower', aspect='auto', cmap='turbo')
plt.colorbar(label='Function Output')
plt.title('Heatmap of Logic Gate 3', fontsize = 15)
plt.xlabel('Y Value', fontsize = 15)
plt.ylabel('X Value', fontsize = 15)
plt.savefig("HeatMap_Logic_Gate_3.png", dpi = 1000)

