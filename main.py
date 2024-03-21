import sys
import os
from matplotlib import pyplot as plt
from termcolor import colored
from TheMostCommonValueEstimate import TheMostCommonValueEstimate
from TheCollisionEstimate import TheCollisionEstimate
from TheMarkovEstimate import TheMarkovEstimate
from TheCompressionEstimate import TheCompressionEstimate
from T_TupleEstimate import T_TupleEstimate
from LRS_Estimate import LRS_Estimate
from MultiMCW_Estimate import MultiMCW_Estimate
from TheMultiMMCPredictionEstimate import TheMultiMMCPredictionEstimate
from TheLagPredictionEstimate import TheLagPredictionEstimate
from TheLZ78YPredictionEstimate import TheLZ78YPredictionEstimate
from bitstring import BitArray

if __name__ == '__main__':
    args = sys.argv

    try:
        dataFile = args[1]
    except:
        print(colored('No data file was provided', 'red'))
        sys.exit(0)

    checkFile = os.path.isfile(dataFile)

    if not checkFile:
        print(colored('Invalid data file was provided', 'red'))
        sys.exit(0)

    try:
        bitsPerSymbol = int([a for a in args if a.isdigit()][0])
    except:
        bitsPerSymbol = 1

    if bitsPerSymbol < 1  or bitsPerSymbol > 8:
        print(colored('Bits per symbol should be between 1 and 8', 'red'))
        sys.exit(0)

    verbose = '-v' in args
    withGraphic = '-g' in args

    with open(dataFile, 'rb') as file:
        bytes = bytearray(file.read())
        strByte = BitArray(bytes)
        #data = [int(b) for b in strByte.bin]
        data = [b & 2 ** bitsPerSymbol - 1 for b in bytes] #take only first bit

        entropies = []

        print(colored('----------- The Most Common Value Estimate -----------', 'blue'))
        entropy = TheMostCommonValueEstimate(data, verbose)
        entropies.append(entropy)
        print('min entropy estimate: ' + str(entropy))

        if bitsPerSymbol == 1:
            
            print(colored('---------------- The Collision Estimate --------------', 'blue'))
            entropy = TheCollisionEstimate(data, verbose)
            entropies.append(entropy)
            print('min entropy estimate: ' + str(entropy))

            print(colored('---------------- The Markov Estimate -----------------', 'blue'))
            entropy = TheMarkovEstimate(data, verbose)
            entropies.append(entropy)
            print('min entropy estimate: ' + str(entropy))
            
            print(colored('-------------- The Compression Estimate --------------', 'blue'))
            entropy = TheCompressionEstimate(data, verbose)
            entropies.append(entropy)
            print('min entropy estimate: ' + str(entropy))
            
        print(colored('------------------ T-Tuple Estimate ------------------', 'blue'))
        entropy = T_TupleEstimate(data, verbose)
        entropies.append(entropy)
        print('min entropy estimate: ' + str(entropy))

        print(colored('-------------------- LRS Estimate --------------------', 'blue'))
        entropy = LRS_Estimate(data, verbose)
        entropies.append(entropy)
        print('min entropy estimate: ' + str(entropy))

        print(colored('----------------- Multi MCW Estimate -----------------', 'blue'))
        entropy = MultiMCW_Estimate(data, verbose)
        entropies.append(entropy)
        print('min entropy estimate: ' + str(entropy))

        print(colored('------------ The Lag Prediction Estimate ------------', 'blue'))
        entropy = TheLagPredictionEstimate(data, verbose)
        entropies.append(entropy)
        print('min entropy estimate: ' + str(entropy))

        print(colored('--------- The Multi MMC Prediction Estimate ---------', 'blue'))
        entropy = TheMultiMMCPredictionEstimate(data, verbose)
        entropies.append(entropy)
        print('min entropy estimate: ' + str(entropy))

        print(colored('----------- The LZ78Y Prediction Estimate -----------', 'blue'))
        entropy = TheLZ78YPredictionEstimate(data, verbose)
        entropies.append(entropy)
        print('min entropy estimate: ' + str(entropy))

        print(colored('result entropy estimate: ' + str(min(entropies)), 'green'))


    labels = [
        'The Most Common Value Estimate',
        'The Collision Estimate',
        'The Markov Estimate',
        'The Compression Estimate',
        'T-Tuple Estimate',
        'LRS Estimate',
        'Multi MCW Estimate',
        'The Lag Prediction Estimate',
        'The Multi MMC Prediction Estimate',
        'The LZ78Y Prediction Estimate',
    ]

    if bitsPerSymbol != 1:
        labels = [
        'The Most Common Value Estimate',
        'T-Tuple Estimate',
        'LRS Estimate',
        'Multi MCW Estimate',
        'The Lag Prediction Estimate',
        'The Multi MMC Prediction Estimate',
        'The LZ78Y Prediction Estimate',
        ]

    colors =  ['#167288','#8cdaec','#b45248','#d48c84','#a89a49', '#d6cfa2', '#3cb464', '#9bddb1', '#643c6a', '#836394']


    if bitsPerSymbol != 1:
        colors =['#167288','#8cdaec','#b45248','#d48c84','#a89a49','#a89a49', '#d6cfa2']
        
    if withGraphic:
        plt.figure(figsize = (len(entropies) + 4, 5))
        plt.bar([i for i in range(len(entropies))], entropies, color = colors)
        handles = [plt.Rectangle((0, 0), 1, 1, color = color) for color in colors]
        plt.legend(handles, labels, loc='upper left', bbox_to_anchor = (1.01, 1), prop = { "size": 12 },)
        plt.tight_layout()
        plt.xticks([])
        plt.show()
        
