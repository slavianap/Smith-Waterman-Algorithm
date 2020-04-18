# -*- coding: utf-8 -*-

"""

LIFE97011 - Computing
Python Programming - Assessed Exercise No. 3
Task: Smith Waterman local alignment
@Author: Slaviana Pavlovich
@CID: 01739756

"""

# Importing Python packages
import numpy as np

# Creating the scoring matrix with the given scores for match, mismatch and gap
def ScoringMatrix(seq1, seq2, MATCH = 1, MISMATCH = -1, GAP = -1):
    # Initialising a NumPy zeros array with the shape of (the first
    # sequence + 1 and the second sequence + 1) and of the integer data type)
    row = len(seq1) + 1
    col = len(seq2) + 1
    M = np.zeros(shape = (row, col), dtype = 'int')
    
    # The score of the cell in the matrix depends on the three cells' values:
    # on the left, above and diagonal (relative to the cell)
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score
            match = MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
            diagonal_score = M[i - 1, j - 1] + match
            # Calculating the score above the cell
            above_score = M[i - 1, j] + GAP
            # Calculating the score on the left
            left_score = M[i, j - 1] + GAP
            # The maximum score is then picked and stored in the cell of the matrix
            M[i, j] = max(0, diagonal_score, above_score, left_score)
    return M
            
# 
#def backtrack()            
    

# Testing    
seq1, seq2 = 'GGTTGACTA', 'TGTTACGG'
L = ScoringMatrix(seq1, seq2)