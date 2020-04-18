# -*- coding: utf-8 -*-

"""

LIFE97011 - Computing
Python Programming - Assessed Exercise No. 3
Task: Smith Waterman (SW) local alignment
@Author: Slaviana Pavlovich
@CID: 01739756

"""

# Importing Python packages
import numpy as np

# Score constants
MATCH = 1
MISMATCH = -1
GAP = -1

# Creating the scoring matrix with the given scores for match, mismatch and gap
def scoring_matrix(seq1, seq2):
    # Initialising a NumPy zeros array with the shape of (the first
    # sequence + 1 and the second sequence + 1) and of the integer data type)
    row = len(seq1) + 1
    col = len(seq2) + 1
    matrix = np.zeros(shape = (row, col), dtype = 'int')    
    
    # The score of the cell in the matrix depends on the three cells' values:
    # on the left, above and diagonal (relative to the cell)
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score
            match_value = MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
            diagonal_score = matrix[i - 1, j - 1] + match_value
            
            # Calculating the score above the cell
            above_score = matrix[i - 1, j] + GAP
            
            # Calculating the score on the left
            left_score = matrix[i, j - 1] + GAP
            
            # The maximum score is then picked and stored in the cell of the matrix
            matrix[i, j] = max(0, diagonal_score, above_score, left_score)
    return matrix
            

#def backtrack()            
  

# Reading the fasta file and keeping the required sequence's name and sequence
def fasta_reader(sequence_file):
    lines = open(sequence_file).readlines()
    sequence_name_row = lines[0][1:]
    sequence = lines[1]
    return sequence_name_row.replace(" ", "").strip(), sequence.strip()
        

# Importing two required fasta sequences
file_1_name, file_1 = fasta_reader("Sequence1.fasta")
file_2_name, file_2  = fasta_reader("Sequence2.fasta")
# Running the Smith Waterman local alignment
M = scoring_matrix(file_1, file_2)