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


def TheMultiMMCPredictionEstimate(S, verbose=False):
    L = len(S)
    k = len(set(S))

    D = 16
    N = L-2
    subpredict = [None for i in range(D)]
    correct = [0 for i in range(N)]
    depths = range(1,D+1)

    M = [{} for i in range(D)]

    scoreboard = [0 for i in range(D)]
    winner = 0

    for i in range(3, L+1):
        if verbose and i%10000 ==0:
            sys.stdout.write("\rComputing MultiMMC Prediction Estimate: %d percents" % (float(i)/L*100))
            sys.stdout.flush()

        for d in depths:
            if d<i-1:
                if M[d-1].get(tuple(S[i-d-2:i-2]),0) == 0:
                    M[d-1][tuple(S[i-d-2:i-2])] = Counter()
                M[d-1][tuple(S[i-d-2:i-2])][S[i-2]] += 1

        for d in depths:
            ymax = 0
            if M[d-1].get(tuple(S[i-d-1:i-1]), None) != None:
                for y in sorted(M[d-1][tuple(S[i-d-1:i-1])].keys()):
                    if M[d-1][tuple(S[i-d-1:i-1])][y] >= M[d-1][tuple(S[i-d-1:i-1])][ymax]:
                        ymax = y
                
                if M[d-1][tuple(S[i-d-1:i-1])][ymax] > 0:
                    subpredict[d-1] = ymax
                    
            else:
                subpredict[d-1] = None

        predict = subpredict[winner]

        if predict == S[i-1]:
            correct[i-2-1] = 1

        for d in range(D):
            if subpredict[d] == S[i-1]:
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