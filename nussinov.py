import sys
import numpy as np
from graphviz import Digraph
#initialize dp matrix
#dp[i][j] = maximal number of bp for subseq from i->j, seek to minimize G

# Background: 
# RNA is made up of {A, U, C, G}
#watson crick {(A, U), (C, G)}
# (G, U) wobble pair for codon-anticodon

MIN_LOOP_LENGTH, SEQUENCE = int(sys.argv[1]), sys.argv[2]
#minimal loop length = minimum unpaired bases in a loop

valid_pairs = set([('A','U'), ('U','A'), ('C', 'G'), ('G','C'), ('G', 'U'), ('U', 'G')])
def isPair(x, y):
    return (x, y) in valid_pairs
#tests
print(isPair('A', 'C'))

def forwardsweep(SEQUENCE):
    N = len(SEQUENCE)
    dp = np.array([[0 for j in range(N)] for i in range(N)])
    #diagonal sweep from major diagonal uprightwards
    for diag in range(2, N+1):
        for i in range(N-diag+1):
            j = i+diag-1
            if(i>j):
                dp[i][j] = 0    
            else:
                #if i or j is unpaired
                IJUnpaired = max(dp[i+1, j], dp[i, j-1])
                # if pairing valid for i<->j
                IJPaired = max(dp[i+1, j], dp[i, j-1], dp[i+1, j-1]+1) if isPair(SEQUENCE[i], SEQUENCE[j]) else 0
                state = lambda i, k: 1 if isPair(i, k) else 0
                #iff pairing  possible for i<->k
                for k in range(i+1, j):
                    Bifurcate = max(dp[i+1, j], dp[i, j-1], dp[i+1, k-1]+dp[k+1, j]+state(i, k))
                    dp[i, j] = max(IJUnpaired, IJPaired, Bifurcate) 

    return dp
print(forwardsweep("GGUCCAC"));

# TRACEBACK
# suboptimals track -delBP 

# Dot notation converter

#connections = []
#def recurse(i, j):
#    if(i>j):
#        return
#    if(dp[i+1, j]==dp[i][j]):
#        connections.append(i, j)
#        recurse(i+1, j)
#    elif(dp[i, j-1]==dp[i][k]):
#        connections.append(i, j)
#        recurse(i+1, j)
