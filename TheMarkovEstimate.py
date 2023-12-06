from math import log2

def TheMarkovEstimate(data, verbose):
    L = len(data)
    bits = 128

    pi = [0, 0]
    P = [[0,0], [0,0]]

    for i in range(L):
        pi[data[i]] += 1

        if(i != L - 1):
            j = i + 1
            P[data[i]][data[j]] += 1

    P[0] = [p_i / sum(P[0]) for p_i in P[0]]
    P[1] = [p_i / sum(P[1]) for p_i in P[1]]

    pi = [pi_i / L for pi_i in pi]

    if(verbose):
        print('pi_0 = ' + str(pi[0]))
        print('pi_1 = ' + str(pi[1]))

    strings = [
        pi[0] * P[0][0] ** 127,
        pi[0] * P[0][1] ** 64 * P[1][0] ** 63,
        pi[0] * P[0][1] * P[1][0] ** 126,
        pi[1] * P[1][0] * P[0][0] ** 126,
        pi[1] * P[1][0] ** 64 * P[0][1] ** 63,
        pi[1] * P[1][1] ** 127,
    ]

    if(verbose):
        print('p_00 = ' + str(P[0][0]))
        print('p_01 = ' + str(P[0][1]))
        print('p_10 = ' + str(P[1][0]))
        print('p_11 = ' + str(P[1][1]))

    p_u = max(strings)

    if(verbose):
        print('p_u = ' + str(p_u))

    H = min(-log2(p_u) / bits, 1)

    return H
    