import math
import numpy


def gaussian(x, mu, sig):
    return (1 / math.sqrt(2 * math.pi * sig ** 2)) * \
        numpy.exp(-(x - mu) ** 2 / (2 * sig ** 2))

def logistic(x, x0, L, M, k):
    return M + (L / (1 + numpy.exp(-k * (x - x0))))
