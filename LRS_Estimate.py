from collections import Counter
from scipy.stats import norm
from math import sqrt, factorial, pow, log

e = 0.995

def find_tuples(s, t):
    return zip(*[s[i:] for i in range(t)])

def find_u(s, threshold):
    count = threshold + 1 
    u = 0
    while count >= threshold:
        u += 1
        c = Counter(find_tuples(s, u))
        count = c.most_common(1)[0][1]
    return u

def find_v(s, u):
    count = 3
    v = u
    while count >= 2:
        v += 1
        c = Counter(find_tuples(s, v))
        count = c.most_common(1)[0][1]
    return v-1

def LRS_Estimate(s, verbose):
    quantile = norm.ppf(e)
    k = 35
    L = len(s)

    u = find_u(s, k)
    if verbose:
        print ("u =",u)

    v = find_v(s, u)
    if verbose:
        print ("v =",v)


    P = []
    for W in range(u, v + 1):
        C = Counter(find_tuples(s, W))
        numerator = 0
        for c in C.values():
            if c == 2:
                numerator += 1
            elif c == 3:
                numerator += 3
            elif c > 3:
                numerator += factorial(c)/2/factorial(c-2)
        
        denom = (L-W+1)*(L-W) / 2 
        P.append(pow(float(numerator)/denom, 1.0/W))
    
    p_max = max(P)
    if verbose:
        print ("p_max:", p_max)    

    p_u = min(1, p_max + quantile * sqrt(p_max * (1 - p_max) / (L - 1)))
    if verbose:
        print ("p_u:", p_u)  

    return -log(p_u,2)
