import sys
import types
import numpy as np
from graphviz import Digraph
#initialize dp matrix
#dp[i][j] = maximal number of bp for subseq from i->j, seek to minimize G

# Background: 
# RNA is made up of {A, U, C, G}
#watson crick {(A, U), (C, G)}
# (G, U) wobble pair for codon-anticodon

#MIN_LOOP_LENGTH, SEQUENCE = int(sys.argv[1]), sys.argv[2]
#minimal loop length = minimum unpaired bases in a loop

valid_pairs = set([('A','U'), ('U','A'), ('C', 'G'), ('G','C'), ('G', 'U'), ('U', 'G')])
def isPair(x, y):
    return (x, y) in valid_pairs
#tests
print(isPair('A', 'C'))


def getDP(dp, i, j, SEQ):
    if(i>j):
        dp[i][j] = 0   
        return 
    else:          
        #if i or j is unpaired
        IJUnpaired = max(dp[i+1, j], dp[i, j-1])
        # if pairing valid for i<->j
        #revert to 0-indexing for the string
        state = lambda x, y: 1 if isPair(SEQ[x-1], SEQ[y-1]) else 0
        IJPaired = dp[i+1, j-1] + state(i, j)
        #iff pairing  possible for i<->k
        if i < j - 1:
            Bifurcate = max([dp[i, k] + dp[k+1, j] for k in range(i, j-1)])
        else: Bifurcate = 0
        dp[i, j] = max(IJUnpaired, IJPaired, Bifurcate) 

def forwardsweep(SEQUENCE):
    N = len(SEQUENCE)
    dp = np.array([[-1 for j in range(N+1)] for i in range(N+1)])
    #diagonal sweep from major diagonal uprightwards
    for i in range(1, N+1):
        dp[i, i] = 0
    for i in range(2, N+1):
        dp[i, i-1] = 0

    for diag in range(2, N+1):
        for j in range(diag, N+1):
            i = j-diag+1
            getDP(dp, i, j, SEQUENCE)
    return dp
print(forwardsweep("GGGAAAUCC"));

# TRACEBACK
#def traceback(dp: ):
#    def recurse(i, j):
#    if(i>j):
#        return
#    if(dp[i+1, j]==dp[i][j]):
#        connections.append(i, j)
#        recurse(i+1, j)
#    elif(dp[i, j-1]==dp[i][k]):
#        connections.append(i, j)
#        recurse(i+1, j)
# suboptimals track -delBP 

# Dot notation converter

connections = []
