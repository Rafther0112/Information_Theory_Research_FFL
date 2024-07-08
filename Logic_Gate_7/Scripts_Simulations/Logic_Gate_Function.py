import numpy as np 

def f(A, B, r, a):
    term1 = r[0] * (a * A) / ((1 + A) * (1 + B))
    term2 = r[1] * (a * B) / ((1 + A) * (1 + B))
    term3 = r[2] * (a * (A * B)) / ((1 + A) * (1 + B))

    result = term1 + term2 - term3

    return result