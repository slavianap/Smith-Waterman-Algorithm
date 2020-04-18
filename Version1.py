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

# Defining the constants for the scores
MATCH = 1
MISMATCH = -1
GAP = -1

# Reading the fasta file and keeping the formatted sequence's name and sequence
def fasta_reader(sequence_file):
    lines = open(sequence_file).readlines()
    sequence_name_row = lines[0][1:]
    sequence = lines[1]
    return sequence_name_row.replace(" ", "").strip(), sequence.strip()

# Creating the scoring matrix
def scoring_matrix(seq1, seq2):
    row = len(seq1) + 1
    col = len(seq2) + 1
    matrix = np.zeros(shape = (row, col), dtype = 'int')   
    
    maximum_score = 0
    maximum_score_location = (0,0)
    
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score (match score)
            match_value = MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
            diagonal_score = matrix[i - 1, j - 1] + match_value
            
            # Calculating the vertical gap score
            vertical_score = matrix[i - 1, j] + GAP
            
            # Calculating the horizontal gap score
            horizontal_score = matrix[i, j - 1] + GAP
            
            # Taking the highest score 
            matrix[i, j] = max(0, diagonal_score, vertical_score, horizontal_score)
        
        # Tracking the cell with the maximum score
        if matrix[i][j] >= maximum_score:
            maximum_score = matrix[i][j]
            maximum_score_location = (i, j) 
            
    return matrix, maximum_score_location


            
#def traceback()            





        
# Importing two required fasta sequences
file_1_name, file_1 = fasta_reader("Sequence1.fasta")
file_2_name, file_2  = fasta_reader("Sequence2.fasta")

# Running the Smith Waterman local alignment
test_matrix, test_maximum_score_location = scoring_matrix(file_1, file_2)