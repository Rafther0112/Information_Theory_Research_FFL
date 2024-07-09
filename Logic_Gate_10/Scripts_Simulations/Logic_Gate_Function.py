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

def f(A, B, r, a):
    term1 = r[0] * (a * A) / ((1 + A) * (1 + B))
    term2 = r[1] * (a * B) / ((1 + A) * (1 + B))
    term3 = r[2] * (a * (A * B)) / ((1 + A) * (1 + B))

    result = a - (term1 + term2 - term3)

    return result