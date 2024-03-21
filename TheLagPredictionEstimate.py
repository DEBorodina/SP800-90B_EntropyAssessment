import sys
import math
from collections import Counter
from decimal import * 
from scipy.stats import norm

e = 0.995

def mostCommon(S,c):
    maxcount = c.most_common()[0][1]
    maxsymb = None

    lastindex = len(S)
    reverse = S[::-1]
    for s in set(S):
        if (c[s] == maxcount) and (reverse.index(s) < lastindex):
            lastindex = reverse.index(s)
            maxsymb = s

    return maxsymb

def calcPavg(C, N, verbose):
    quantile = norm.ppf(e)

    p_global = float(C)/N

    if(p_global ==0 ):
        p_globalprime = 1 - 0.01 ** (1 / N)
    else:
        p_globalprime = p_global + quantile * math.sqrt(float(p_global) * (1 - p_global) / (N - 1))

    if verbose:
        print("\nP_global: %f" % p_globalprime)

    return p_globalprime

def calc_qn(p,r,n):
    p = Decimal(str(p))
    q = 1-p
    x = Decimal(str(find_root(p,r)))
    
    qn = (1-p*x)/Decimal((r+1-r*x)*q)
    qn = qn/(x**(n+1))

    return qn

def find_root(p,r):
    p = Decimal(str(p))
    q = 1-p

    s = Decimal(1)
    for i in range(10):
            s = 1+q*(p**r)*(s**(r+1))
    return s

def findMaxRun(correct):
    run = 0
    maxrun = 0
    for i in correct:
        if i == 0:
            if run > maxrun:
                maxrun = run
            run = 0
        elif i==1:
            run += 1
        else:
            raise ValueError("correct array contains non-binary values")
    if run > maxrun:
        maxrun = run

    return maxrun

def calcRun(correct):
    N = len(correct)
    alpha = 0.99
       
    r = findMaxRun(correct)
    alpha = Decimal(str(alpha))
    
    p = Decimal(str(0.5))
    adj = Decimal(str(0.5))

    qn = calc_qn(p,r+1,N)
    
    for _ in range(30): 
            adj /= 2
            if qn > alpha:
                    p += adj
            else:
                    p -= adj
                    
            qn = calc_qn(p,r+1,N)
            if abs(qn-alpha) <= 0.0001: break
   
    return p

def TheLagPredictionEstimate(S, verbose=False):
    L = len(S)
    k = len(set(S))

    D = 128
    N = L - 1
    correct = [0 for _ in range(N)]
    lag = [None for _ in range(D)]

    scoreboard = [0 for _ in range(D)]
    winner = 0

    for i in range(2, L + 1):
        if verbose and i % 1000 == 0:
            sys.stdout.write("\rComputing Lag Prediction Estimate: %d percents" % (float(i)/L*100))
            sys.stdout.flush()

        for d in range(1,D + 1):
            if d < i:
                lag[d - 1] = S[i - d - 1]
            else:
                lag[d - 1] = None

        prediction = lag[winner]

        if prediction == S[i - 1]:
            correct[i - 2] = 1

        for d in range(D):
            if lag[d] == S[i - 1]:
                scoreboard[d] += 1
                if scoreboard[d] >= scoreboard[winner]:
                    winner = d
    
    C = sum(correct)

    P_global = calcPavg(C, N, verbose)

    P_local= calcRun(correct)

    minH = -math.log(max(P_global, P_local, 1/k),2)
    
    if verbose:
        print("P_global': %f" % P_global)
        print("P_local: %f"% P_local)

    return minH