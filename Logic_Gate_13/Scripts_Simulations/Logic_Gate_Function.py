import numpy as np

def f(R, A, K, a):
    # Ensure K is a numpy array for element-wise operations
    K = np.array(K)

    term1 = K[3] * (K[1] + R * K[4])
    term2 = K[3] * (A * K[0] + K[1] + R * K[4])
    term3 = (R * K[2] + K[3]) * (K[1] + R * K[4])
    term4 = A * K[0] * (K[3] + R * K[4])

    numerator = term1 - term2
    denominator = term3 + term4
    fraction = numerator / denominator

    result = a * (1 + fraction)

    return result