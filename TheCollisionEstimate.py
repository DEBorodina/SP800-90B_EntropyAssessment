from statistics import mean, stdev
from scipy.stats import norm
from math import sqrt, log2, exp
import mpmath

e = 0.995 

def F(q):
    k = 3
    z = 1 / q

    return mpmath.gammainc(k,z) * z ** (-k) * exp(z)

def function(p):
    q = 1 - p

    firstAdd = p * q ** (-2) * (1 + (p ** (-1) - q ** (-1)) /  2)
    secondAdd = p * q ** (-1) * (p ** (-1) - q ** (-1)) / 2

    return firstAdd * F(q) - secondAdd

def binarySearch(u_hat):
    l = 0.5
    r = 1
    step = (r - l)
    p = step / 2 + l

    eps = 10 ** (-8)

    if u_hat > function(l):
        return False, 0

    funVal = function(p)
    while abs(u_hat - funVal) > eps:
        step /= 2
        if u_hat < funVal:
            p += step
            if p == 1:
                p -= 0.0001
        else:
            p -= step
            if p < l:
                p = l
        funVal = function(p)

    return True, p

def TheCollisionEstimate(data, verbose):
    quantile = norm.ppf(e)

    t = []

    index = 0
    for j, val in enumerate(data):
        if val in data[index:j]:
            t.append(j + 1 - index)
            index = j + 1

    v = len(t)

    if verbose:
        print('v = ' + str(v))

    mu = mean(t)
    s = stdev(t)

    if verbose:
        print('mu = ' + str(mu))
        print('s = ' + str(s))

    u_hat = mu - quantile * s / sqrt(v)

    if verbose:
        print('u_hat = ' + str(u_hat))

    valid, p = binarySearch(u_hat)

    if verbose:
        print('p = ' + str(p))
    
    if not valid:
        H = 1
    else:
        H = -log2(p)

    return H 


