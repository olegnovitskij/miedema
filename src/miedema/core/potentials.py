import numpy as np


def Morse_potential(x, n=12, m=6):
    gamma = (n * m / 2) ** 0.5
    return np.exp(-2 * gamma * (x - 1)) - 2 * np.exp(-gamma * (x - 1))


def L_J_potential(x, n=12, m=6):
    return (m / (n - m)) * (x ** (-n) - (n / m) * x ** (-m))
