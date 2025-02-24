#initialize dp matrix
#dp[i][j] = maximal number of bp for subseq from i->j

# RNA is made up of {A, U, C, G}
#watson crick {(A, U), (C, G)}
# (G, U) wobble pair

#minimal loop length = minimum unpaired bases in a loop

#dp[i, j]= 
# 0 if i>j
# max(dp[i+1, j], dp[i, j-1]) and at least one unpaired
# max(dp[i+1, j], dp[i, j-1], dp[i+1, j+1]+1), pair i<->j
# for all j, i<k<j
#   state = 1 if i and k pair else 0
#   max(dp[i+1, k-1]+dp[k+1, j]) if i pairs with k