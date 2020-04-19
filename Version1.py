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

# Implementing the Smith Waterman (SW) local alignment
def SW(seq1, seq2):
    maximum_score = 0
    STOP, LEFT, UP, DIAGONAL = range(4)
    aligned_seq1 = ""
    aligned_seq2 = ""    
    
    # Creating the empty matrices for storing scores and tracing
    row = len(seq1) + 1
    col = len(seq2) + 1
    matrix = np.zeros(shape = (row, col), dtype = 'int')  
    tracing_matrix = np.zeros(shape = (row, col), dtype = 'int')  
    
    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score (match score)
            match_value = MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
            diagonal_score = matrix[i - 1, j - 1] + match_value
            # Calculating the vertical gap score
            vertical_score = matrix[i - 1, j] + GAP
            # Calculating the horizontal gap score
            horizontal_score = matrix[i, j - 1] + GAP
            
            # Tracking where the cell's value is coming from    
            if matrix[i][j] == 0: 
                tracing_matrix[i][j] = STOP
            elif matrix[i][j] == horizontal_score: 
                tracing_matrix[i][j] = LEFT
            elif matrix[i][j] == vertical_score: 
                tracing_matrix[i][j] = UP
            else: 
                tracing_matrix[i][j] = DIAGONAL            
            
            # Taking the highest score 
            matrix[i, j] = max(0, diagonal_score, vertical_score, horizontal_score)
         
        # Tracking the cell with the maximum score
        if matrix[i][j] >= maximum_score:
            maximum_i = i
            maximum_j = j
            maximum_score = matrix[i][j]
           
    # Assigning to i and j the cell's position with the highest score
    i, j = maximum_i, maximum_j
    
    while tracing_matrix[i, j] != STOP:
        if tracing_matrix[i, j] == DIAGONAL:
            current_aligned_seq1 = seq1[i - 1]
            current_aligned_seq2 = seq2[j - 1]
            i = i - 1
            j = j - 1
        elif tracing_matrix[i, j] == UP:
            current_aligned_seq1 = seq1[i - 1]
            current_aligned_seq2 = '-'
            i = i - 1            
        else:
            current_aligned_seq1 = '-'
            current_aligned_seq2 = seq2[j - 1]
            j = j - 1
        aligned_seq1 = aligned_seq1 + current_aligned_seq1
        aligned_seq2 = aligned_seq2 + current_aligned_seq2
    
    # Reversing the order of the sequences
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]
    
    return(aligned_seq1, aligned_seq2)

# Importing two required fasta sequences
file_1_name, file_1 = fasta_reader("Sequence1.fasta")
file_2_name, file_2  = fasta_reader("Sequence2.fasta")

# Running the Smith Waterman local alignment and displaying the alignment
output_1, output_2 = SW(file_1, file_2)
print(file_1_name + ' ' + output_1 + '\n' + file_2_name + ' ' + output_2)