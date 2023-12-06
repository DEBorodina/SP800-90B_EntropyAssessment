import sys
from matplotlib import pyplot as plt
from termcolor import colored
from TheMostCommonValueEstimate import TheMostCommonValueEstimate
from TheCollisionEstimate import TheCollisionEstimate
from TheMarkovEstimate import TheMarkovEstimate
from TheCompressionEstimate import TheCompressionEstimate
from T_TupleEstimate import T_TupleEstimate

if __name__ == '__main__':
    args = sys.argv

    try:
        dataFile = args[1]
    except:
        print(colored('No data file was provided', 'red'))
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

        print(colored('result entropy estimate: ' + str(min(entropies)), 'green'))

    labels = [
        'The Most Common Value Estimate',
        'The Collision Estimate',
        'The Markov Estimate',
        'The Compression Estimate',
        'T-Tuple Estimate',
    ]

    if bitsPerSymbol != 1:
        labels = [
        'The Most Common Value Estimate',
        'T-Tuple Estimate',
        ]

    colors = ['#845EC2','#D65DB1','#FF9671','#FFC75F','#F9F871']

    if withGraphic:
        plt.figure(figsize = (len(entropies) + 4, 5))
        plt.bar([i for i in range(len(entropies))], entropies, color = colors)
        handles = [plt.Rectangle((0, 0), 1, 1, color = color) for color in colors]
        plt.legend(handles, labels, loc='upper left', bbox_to_anchor = (1.01, 1))
        plt.tight_layout()
        plt.xticks([])
        plt.show()
        
