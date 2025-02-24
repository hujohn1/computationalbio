import sys
from graphviz import Digraph
#initialize dp matrix
#dp[i][j] = maximal number of bp for subseq from i->j, seek to minimize G
#very similar to landau levenshtein mods

# Background: 
# RNA is made up of {A, U, C, G}
#watson crick {(A, U), (C, G)}
# (G, U) wobble pair for codon-anticodon

MIN_LOOP_LENGTH, SEQUENCE = sys.argv[0], sys.argv[1]
#minimal loop length = minimum unpaired bases in a loop

N = len(SEQUENCE)
dp = [[] * N]
[dp[i][j] == 0 if i==j else 0 for i in range (0, N) for j in range(0, N)]


valid_pairs = set([('A','U'), ('U','A'), ('C', 'G'), ('G','C'), ('G', 'U'), ('U', 'G')])
def isPair(x, y):
    return (x, y) in valid_pairs

i = 1
j = 10
if i>j:
    dp[i][j] = 0
    
else:
    #if i or j is unpaired
    IJUnpaired = dp[i+1, j], dp[i, j-1]
    # if pairing valid for i<->j
    IJPaired = max(dp[i+1, j], dp[i, j-1], dp[i+1, j-1]+1)
    state = lambda i, k: 1 if isPair(i, k) else 0
    #iff pairing  possible for i<->k
    for k in range(i+1, j):
        Bifurcate = max(dp[i+1, j], dp[i, j-1], dp[i+1, k-1]+dp[k+1, j]+state(i, k))
    dp[i, j] = max(IJUnpaired, IJPaired, Bifurcate) 

# TRACEBACK
# suboptimals track -delBP 