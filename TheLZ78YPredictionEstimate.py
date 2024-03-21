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

def TheLZ78YPredictionEstimate(S, verbose=False):
    L = len(S)
    l = len(set(S))

    B = 16
    N = L-B-1
    correct = [0 for i in range(N)]
    maxDictionarySize = 65536

    D = dict()
    dictionarySize = 0

    for i in range(B + 2, L + 1):

        if verbose and i % 10000==0:
            sys.stdout.write("\rComputing LZ78Y Prediction Estimate: %d percents" % (float(i)/L*100))
            sys.stdout.flush()

        for j in range(B, 0, -1):
            k = tuple(S[i - j - 2 :i - 2]) 
            if k not in D and dictionarySize < maxDictionarySize:
                D[k] = dict()
                dictionarySize = dictionarySize + 1
            if k in D: 
                D[k][S[i - 2]] = D[k].get(S[i - 2], 0) + 1

        maxcount = 0
        predict = None
        for j in range(B, 0, -1):
            prev = tuple(S[i - j - 1:i - 1])
            if prev in D:
                for y in sorted(D[prev].keys(), reverse = True):
                    if D[prev][y] > maxcount:
                        predict = y
                        maxcount = D[prev][y]
                
        if predict == S[i - 1]:
            correct[i - B - 1 - 1] += 1
        
    C = sum(correct)

    P_global = calcPavg(C, N, verbose)

    P_local= calcRun(correct)

    minH = -math.log(max(P_global, P_local, 1/l),2)
    
    if verbose:
        print("P_global': %f" % P_global)
        print("P_local: %f"% P_local)

    return minH