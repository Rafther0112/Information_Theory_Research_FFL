import numpy as np
from tqdm import tqdm
from numba import jit,njit
import pandas as pd
import json

@njit
def Hill_Activation(cantidad, sensitivity, expresion_level, hill):
    return expresion_level*((cantidad**hill)/(cantidad**hill + sensitivity**hill))

@njit
def Hill_Represion(cantidad, sensitivity, expresion_level, hill):
    return expresion_level*((sensitivity**hill)/(cantidad**hill + sensitivity**hill))

@njit
def logic_gate_function_13(R, A, K):
    # Ensure K is a numpy array for element-wise operations
    K = np.array(K)

    term1 = K[3] * (K[1] + R * K[4])
    term2 = K[3] * (A * K[0] + K[1] + R * K[4])
    term3 = (R * K[2] + K[3]) * (K[1] + R * K[4])
    term4 = A * K[0] * (K[3] + R * K[4])

    numerator = term1 - term2
    denominator = term3 + term4
    fraction = numerator / denominator

    result = 1 + fraction

    return result