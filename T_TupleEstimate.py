from scipy.stats import norm
from math import sqrt, log2
from collections import Counter

e = 0.995

def T_TupleEstimate(data, verbose):
    quantile = norm.ppf(e)
    k = 35 #3

    n = len(data)

    Q = []
    i = 0

    while True:
        c = Counter(zip(*[data[j:] for j in range(i + 1)]))
        mostCommonCount = c.most_common(1)[0][1]

        if mostCommonCount < k:
            break

        Q.append(mostCommonCount)
        i += 1
        
    t = i

    if(verbose):
        print('t = ' + str(t))
    
    P = []
    for i in range(t):
        P.append(Q[i] / (n - (i + 1) + 1))
        P[i] = P[i] ** (1 / (i + 1))
    
    p_max = max(P)
    p_u = min(1, p_max + quantile * sqrt(p_max * (1 - p_max) / (n - 1)))

    if(verbose):
        print('p_max = ' + str(p_max))
        print('p_u = ' + str(p_u))
    
    H = -log2(p_u)

    return H

