from collections import Counter
from scipy.stats import norm
from math import sqrt, log2

e = 0.995

def TheMostCommonValueEstimate(data, verbose):
    L = len(data)
    quantile = norm.ppf(e)

    _, occurrence_count = Counter(data).most_common(1)[0]
    p_hat = occurrence_count / L

    if(verbose):
        print('p = ' + str(p_hat))

    p_u = min(1, p_hat + quantile * sqrt(p_hat * (1 - p_hat) / (L - 1)))

    if(verbose):
        print('p_u = ' + str(p_u))

    H = -log2(p_u)

    return H


    
