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
def logic_gate_function_10(B, A, K, a):
    term1 = K[0] *  A / ((K[3] + A) * (K[4] + B))
    term2 = K[1] * B / ((K[3] + A) * (K[4] + B))
    term3 = (K[2] * A * B) / ((K[3] + A) * (K[4] + B))

    result = a*(1 - (term1 + term2 - term3))

    return result