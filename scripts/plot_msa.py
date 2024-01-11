import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
import sys
import os

# Check that two command-line arguments have been provided
if len(sys.argv) != 2:
    print("USAGE: python plot_msa.py <input_msa.a3m>", file=sys.stderr)
    print("OUTPUT:")
    print("   input_msa_msa.png     : png file with the same name as the input")
    sys.exit(1)

def blosum62_score(res1, res2):
    blosum62 = { 
    ('A', 'A'): 4, ('R', 'A'):  -1, ('N', 'A'):  -2, ('D', 'A'):  -2, ('C', 'A'):  0, ('Q', 'A'):  -1, ('E', 'A'):  -1, ('G', 'A'):  0, ('H', 'A'):  -2, ('I', 'A'):  -1, ('L', 'A'):  -1, ('K', 'A'):  -1, ('M', 'A'):  -1, ('F', 'A'):  -2, ('P', 'A'):  -1, ('S', 'A'):  1, ('T', 'A'):  0, ('W', 'A'):  -3, ('Y', 'A'):  -2, ('V', 'A'):  0, 
    ('A', 'R'):  -1, ('R', 'R'):  5, ('N', 'R'):  0, ('D', 'R'):  -2, ('C', 'R'):  -3, ('Q', 'R'):  1, ('E', 'R'):  0, ('G', 'R'):  -2, ('H', 'R'):  0, ('I', 'R'):  -3, ('L', 'R'):  -2, ('K', 'R'):  2, ('M', 'R'):  -1, ('F', 'R'):  -3, ('P', 'R'):  -2, ('S', 'R'):  -1, ('T', 'R'):  -1, ('W', 'R'):  -3, ('Y', 'R'):  -2, ('V', 'R'):  -3, 
    ('A', 'N'):  -2, ('R', 'N'):  0, ('N', 'N'):  6, ('D', 'N'):  1, ('C', 'N'):  -3, ('Q', 'N'):  0, ('E', 'N'):  0, ('G', 'N'):  0, ('H', 'N'):  1, ('I', 'N'):  -3, ('L', 'N'):  -3, ('K', 'N'):  0, ('M', 'N'):  -2, ('F', 'N'):  -3, ('P', 'N'):  -2, ('S', 'N'):  1, ('T', 'N'):  0, ('W', 'N'):  -4, ('Y', 'N'):  -2, ('V', 'N'):  -3, 
    ('A', 'D'):  -2, ('R', 'D'):  -2, ('N', 'D'):  1, ('D', 'D'):  6, ('C', 'D'):  -3, ('Q', 'D'):  0, ('E', 'D'):  2, ('G', 'D'):  -1, ('H', 'D'):  -1, ('I', 'D'):  -3, ('L', 'D'):  -4, ('K', 'D'):  -1, ('M', 'D'):  -3, ('F', 'D'): -3, ('P', 'D'): -1, ('S', 'D'): 0, ('T', 'D'): -1, ('W', 'D'): -4, ('Y', 'D'): -3, ('V', 'D'): -3,
    ('A', 'C'): 0, ('R', 'C'): -3, ('N', 'C'): -3, ('D', 'C'): -3, ('C', 'C'): 9, ('Q', 'C'): -3, ('E', 'C'): -4, ('G', 'C'): -3, ('H', 'C'): -3, ('I', 'C'): -1, ('L', 'C'): -1, ('K', 'C'): -3, ('M', 'C'): -1, ('F', 'C'): -2, ('P', 'C'): -3, ('S', 'C'): -1, ('T', 'C'): -1, ('W', 'C'): -2, ('Y', 'C'): -2, ('V', 'C'): -1,
    ('A', 'Q'): -1, ('R', 'Q'): 1, ('N', 'Q'): 0, ('D', 'Q'): 0, ('C', 'Q'): -3, ('Q', 'Q'): 5, ('E', 'Q'): 2, ('G', 'Q'): -2, ('H', 'Q'): 0, ('I', 'Q'): -3, ('L', 'Q'): -2, ('K', 'Q'): 1, ('M', 'Q'): 0, ('F', 'Q'): -3, ('P', 'Q'): -1, ('S', 'Q'): 0, ('T', 'Q'): -1, ('W', 'Q'): -2, ('Y', 'Q'): -1, ('V', 'Q'): -2,
    ('A', 'E'): -1, ('R', 'E'): 0, ('N', 'E'): 0, ('D', 'E'): 2, ('C', 'E'): -4, ('Q', 'E'): 2, ('E', 'E'): 5, ('G', 'E'): -2, ('H', 'E'): 0, ('I', 'E'): -3, ('L', 'E'): -3, ('K', 'E'): 1, ('M', 'E'): -2, ('F', 'E'): -3, ('P', 'E'): -1, ('S', 'E'): 0, ('T', 'E'): -1, ('W', 'E'): -3, ('Y', 'E'): -2, ('V', 'E'): -2,
    ('A', 'G'): 0, ('R', 'G'): -2, ('N', 'G'): 0, ('D', 'G'): -1, ('C', 'G'): -3, ('Q', 'G'): -2, ('E', 'G'): -2, ('G', 'G'): 6, ('H', 'G'): -2, ('I', 'G'): -4, ('L', 'G'): -4, ('K', 'G'): -2, ('M', 'G'): -3, ('F', 'G'): -3, ('P', 'G'): -2, ('S', 'G'): 0, ('T', 'G'): -2, ('W', 'G'): -2, ('Y', 'G'): -3, ('V', 'G'): -3,
    ('A', 'H'): -2, ('R', 'H'): 0, ('N', 'H'): 1, ('D', 'H'): -1, ('C', 'H'): -3, ('Q', 'H'): 0, ('E', 'H'): 0, ('G', 'H'): -2, ('H', 'H'): 8, ('I', 'H'): -3, ('L', 'H'): -3, ('K', 'H'): -1, ('M', 'H'): -2, ('F', 'H'): -1, ('P', 'H'): -2, ('S', 'H'): -1, ('T', 'H'): -2, ('W', 'H'): -2, ('Y', 'H'): 2, ('V', 'H'): -3,
    ('A', 'I'): -1, ('R', 'I'): -3, ('N', 'I'): -3, ('D', 'I'): -3, ('C', 'I'): -1, ('Q', 'I'): -3, ('E', 'I'): -3, ('G', 'I'): -4, ('H', 'I'): -3, ('I', 'I'): 4, ('L', 'I'): 2, ('K', 'I'): -3, ('M', 'I'): 1, ('F', 'I'): 0, ('P', 'I'): -3, ('S', 'I'): -2, ('T', 'I'): -1, ('W', 'I'): -3, ('Y', 'I'): -1, ('V', 'I'): 3,
    ('A', 'L'): -1, ('R', 'L'): -2, ('N', 'L'): -4, ('D', 'L'): -4, ('C', 'L'): -1, ('Q', 'L'): -2, ('E', 'L'): -3, ('G', 'L'): -4, ('H', 'L'): -3, ('I', 'L'): 2, ('L', 'L'): 4, ('K', 'L'): -2, ('M', 'L'): 2, ('F', 'L'): 0, ('P', 'L'): -3, ('S', 'L'): -2, ('T', 'L'): -1, ('W', 'L'): -2, ('Y', 'L'): -1, ('V', 'L'): 1,
    ('A', 'K'): -1, ('R', 'K'): 2, ('N', 'K'): 0, ('D', 'K'): 0, ('C', 'K'): -3, ('Q', 'K'): 1, ('E', 'K'): 1, ('G', 'K'): -2, ('H', 'K'): -1, ('I', 'K'): -3, ('L', 'K'): -2, ('K', 'K'): 5, ('M', 'K'): -1, ('F', 'K'): -3, ('P', 'K'): -1, ('S', 'K'): 0, ('T', 'K'): -1, ('W', 'K'): -3, ('Y', 'K'): -2, ('V', 'K'): -2,
    ('A', 'M'): -1, ('R', 'M'): -1, ('N', 'M'): -2, ('D', 'M'): -3, ('C', 'M'): -1, ('Q', 'M'): 0, ('E', 'M'): -2, ('G', 'M'): -3, ('H', 'M'): -2, ('I', 'M'): 1, ('L', 'M'): 2, ('K', 'M'): -1, ('M', 'M'): 5, ('F', 'M'): 0, ('P', 'M'): -2, ('S', 'M'): -1, ('T', 'M'): -1, ('W', 'M'): -1, ('Y', 'M'): -1, ('V', 'M'): 1,
    ('A', 'F'): -2, ('R', 'F'): -3, ('N', 'F'): -3, ('D', 'F'): -3, ('C', 'F'): -2, ('Q', 'F'): -3, ('E', 'F'): -3, ('G', 'F'): -3, ('H', 'F'): -1, ('I', 'F'): 0, ('L', 'F'): 0, ('K', 'F'): -3, ('M', 'F'): 0, ('F', 'F'): 6, ('P', 'F'): -4, ('S', 'F'): -2, ('T', 'F'): -2, ('W', 'F'): 1, ('Y', 'F'): 3, ('V', 'F'): -1,
    ('A', 'P'): -1, ('R', 'P'): -2, ('N', 'P'): -2, ('D', 'P'): -1, ('C', 'P'): -3, ('Q', 'P'): -1, ('E', 'P'): -1, ('G', 'P'): -2, ('H', 'P'): -2, ('I', 'P'): -3, ('L', 'P'): -3, ('K', 'P'): -1, ('M', 'P'): -2, ('F', 'P'): -4, ('P', 'P'): 7, ('S', 'P'): -1, ('T', 'P'): -1, ('W', 'P'): -4, ('Y', 'P'): -3, ('V', 'P'): -2,
    ('A', 'S'): 1, ('R', 'S'): -1, ('N', 'S'): 1, ('D', 'S'): 0, ('C', 'S'): -1, ('Q', 'S'): 0, ('E', 'S'): 0, ('G', 'S'): 0, ('H', 'S'): -1, ('I', 'S'): -2, ('L', 'S'): -2, ('K', 'S'): 0, ('M', 'S'): -1, ('F', 'S'): -2, ('P', 'S'): -1, ('S', 'S'): 4, ('T', 'S'): 1, ('W', 'S'): -3, ('Y', 'S'): -2, ('V', 'S'): -1,
    ('A', 'T'): 0, ('R', 'T'): -1, ('N', 'T'): 0, ('D', 'T'): -1, ('C', 'T'): -1, ('Q', 'T'): -1, ('E', 'T'): -1, ('G', 'T'): -2, ('H', 'T'): -2, ('I', 'T'): -1, ('L', 'T'): -1, ('K', 'T'): -1, ('M', 'T'): -1, ('F', 'T'): -2, ('P', 'T'): -1, ('S', 'T'): 1, ('T', 'T'): 5, ('W', 'T'): -2, ('Y', 'T'): -2, ('V', 'T'): 0,
    ('A', 'W'): -3, ('R', 'W'): -3, ('N', 'W'): -4, ('D', 'W'): -4, ('C', 'W'): -2, ('Q', 'W'): -2, ('E', 'W'): -3, ('G', 'W'): -2, ('H', 'W'): -2, ('I', 'W'): -3, ('L', 'W'): -2, ('K', 'W'): -3, ('M', 'W'): -1, ('F', 'W'): 1, ('P', 'W'): -4, ('S', 'W'): -3, ('T', 'W'): -2, ('W', 'W'): 11, ('Y', 'W'): 2, ('V', 'W'): -3,
    ('A', 'Y'): -2, ('R', 'Y'): -2, ('N', 'Y'): -2, ('D', 'Y'): -3, ('C', 'Y'): -2, ('Q', 'Y'): -1, ('E', 'Y'): -2, ('G', 'Y'): -3, ('H', 'Y'): 2, ('I', 'Y'): -1, ('L', 'Y'): -1, ('K', 'Y'): -2, ('M', 'Y'): -1, ('F', 'Y'): 3, ('P', 'Y'): -3, ('S', 'Y'): -2, ('T', 'Y'): -2, ('W', 'Y'): 2, ('Y', 'Y'): 7, ('V', 'Y'): -1,
    ('A', 'V'): 0, ('R', 'V'): -3, ('N', 'V'): -3, ('D', 'V'): -3, ('C', 'V'): -1, ('Q', 'V'): -2, ('E', 'V'): -2, ('G', 'V'): -3, ('H', 'V'): -3, ('I', 'V'): 3, ('L', 'V'): 1, ('K', 'V'): -2, ('M', 'V'): 1, ('F', 'V'): -1, ('P', 'V'): -2, ('S', 'V'): -1, ('T', 'V'): 0, ('W', 'V'): -3, ('Y', 'V'): -1, ('V', 'V'): 4
    }
    try:
        return blosum62[(res1, res2)]
    except KeyError:
        if res2 == "-":
            return np.nan
        return 0

def plot_msa(msa_file):
    
    with open(msa_file, "r") as file_read:
        for i, line in enumerate(file_read):
            if i == 0:
                cardinality = line.lstrip("#").replace("\t", "   ")
            else:
                break
    
    header_to_plot = msa_file.split("/")[-1].replace(".a3m", "").replace("__vs__", ":") + "\n"
    
    # Remove lowercase letters from file and save it as temporal
    with open(msa_file, 'r') as f:
        text = f.read()
    new_text = ''.join(c for c in text if not c.islower())
    new_filename = msa_file + '_temporal'
    with open(new_filename, 'w') as f:
        f.write(new_text)

    # Read in the MSA file using Biopython
    alignment = AlignIO.read(new_filename, "fasta")
    
    # Get the number of sequences and the length of the alignment
    num_seqs = len(alignment)
    alignment_len = alignment.get_alignment_length()
    
    # Create a matrix to hold the similarity scores
    similarity_matrix = np.zeros((num_seqs, alignment_len))
    
    # Calculate the similarity scores for each position
    for i in range(alignment_len):
        query_res = alignment[0][i]
        for j in range(num_seqs):
            other_res = alignment[j][i]
            similarity_matrix[j,i] = blosum62_score(query_res, other_res)
    
    # Plot the similarity matrix as a heatmap
    plt.imshow(similarity_matrix, cmap='rainbow', vmin=-4, vmax=9)
    plt.gca().set_aspect('auto')
    plt.colorbar(label='BLOSUM62 score to query position')
    
    # Labels
    plt.title(header_to_plot)
    plt.suptitle(cardinality, y=0.92, x = 0.44, fontsize = 9)
    plt.xlabel('Position')
    plt.ylabel('Sequences')
    
    # Save the plot
    out_file = str(msa_file).split("/")[-1]
    out_file = out_file.replace(".a3m", "_msa.png")
    plt.savefig(out_file, dpi=300)

    # Remove temporal file	
    os.remove(new_filename)

a3m_file=sys.argv[1]
plot_msa(a3m_file)
