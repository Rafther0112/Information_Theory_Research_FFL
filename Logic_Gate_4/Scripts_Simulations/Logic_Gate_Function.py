import numpy as np

def f(A, R, K):
    # Ensure K is a numpy array for element-wise operations
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