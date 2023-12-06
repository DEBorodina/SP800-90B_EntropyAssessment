from math import floor, log2, sqrt
from scipy.stats import norm

def G(z, v, d):
    n = v + d

    s = 0
    common = z ** 2 * sum([log2(u) * (1 - z) **(u - 1) for u in range(1,(d + 1))])
    s += common * v

    newAddOnEveryStep = [log2(u) * (1 - z) **(u - 1)  for u in range((d + 1), n + 1)]

    s += z ** 2 * sum([(n - i - (d + 1)) * newAddOnEveryStep[i] for i in range(n - (d + 1))])

    s += z * sum(newAddOnEveryStep)

    return s / v

def function(p, v, d, b):
    q = (1 - p)/(2 ** b - 1)
    
    return G(p, v, d) + (2 ** b - 1) * G(q, v, d)

def binarySearch(mu_hat, v, d, b):
    l = 2 ** (-b)
    r = 1
    step = (r - l)
    p = step / 2 + l

    eps = 10 ** (-6)

    if mu_hat > function(l, v, d, b):
        return False, 0

    funVal = function(p, v, d, b)

    while abs(mu_hat - funVal) > eps:
        step /= 2.0
        if mu_hat > funVal:
            p -= step
        else:
            p += step
        funVal = function(p, v, d, b)

    return True, p

def TheCompressionEstimate(data, verbose):
    n = len(data)
    b = 6
    n1 = floor(n / b)

    x = [data[i:i + b] for i in range(0, n1 * b, b)]

    d = 1000 #4
    v = n1 - d

    if(verbose):
        print('v = ' + str(v))

    dict = {}
    for i in range(d):
        dict[str(x[i])] = i + 1

    D = [ 0 for _ in range(v)]
    for i in range(d, n1):
        if str(x[i]) in dict:
            D[i - d] = i + 1 - dict[str(x[i])]
        else:
            D[i - d] = i + 1

        dict[str(x[i])] = i + 1

    c = 0.5907
    quantile = norm.ppf(0.995)

    mu = sum([log2(D_i) for D_i in D]) / v
    s = c * sqrt((1 / (v - 1)) * sum([log2(D_i) ** 2 for D_i in D]) - mu ** 2)
    mu_hat = mu - quantile * s / sqrt(v)

    if(verbose):
        print('mu = ' + str(mu))
        print('s = ' + str(s))
        print('mu_hat = ' + str(mu_hat))

    valid, p = binarySearch(mu_hat, v, d, b)

    if(verbose):
        print('p = ' + str(p))

    if not valid:
        H = 1
    else:
        H = -log2(p) / b

    return H 
