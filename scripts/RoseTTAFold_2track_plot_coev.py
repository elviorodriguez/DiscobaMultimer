# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Check if command-line arguments have been provided
if len(sys.argv) != 5:
    print("ERROR: missing positional arguments", file=sys.stderr)
    print("USAGE: python RoseTTAFold_2track_plot_coev.py <npz_file_d> <npz_file_r> <top> <a3m_seqs_N>", file=sys.stderr)
    print("")
    print("   npz_file_d  .npz file with coevolutionary info produced with")
    print("               RoseTTAFold_2track_run.sh. It is expected to have")
    print("               the format <proteinID1__vs__proteinID2.npz>")
    print("   npz_file_r  Reversed .npz file coevolutionary info produced with")
    print("               RoseTTAFold_2track_run.sh. It is expected to have")
    print("               the format <proteinID2__vs__proteinID1.npz>")
    print("   top         Number of contacts to display (integer)")
    print("   a3m_seqs_N  Number of sequences in the original MSA")    
    print("")
    print("OUTPUT:")
    print("")
    print("    1) method-proteinID1__vs__proteinID2-coevolution.png:")
    print("        predicted contacts plot as a heatmap")
    print("        method prefix correspond to the normalization type (direct,")
    print("        reversed, min, max, mean)")
    print("")
    print("    2) method-proteinID1__vs__proteinID2.contacts:")
    print("        tsv file with residues predicted to be in contact and their")
    print("        contact probability")
    print("")
    print("    3) method-proteinID1__vs__proteinID2.chimeraX:")
    print("        commands for ChimeraX to display predicted contacs in the")
    print("        structural models")
    print("")
    print("    4) proteinID1__vs__proteinID2.metrics")
    print("        tsv file with the computed metrics to predict PPIs")
    sys.exit(1)

# Assign the command-line arguments to variables
npz_file_d = sys.argv[1]    # direct
npz_file_r = sys.argv[2]    # reversed
top = int(sys.argv[3])
a3m_seqs_N = int(sys.argv[4])

# Useful information and names
npz_file_dirname = os.path.dirname(npz_file_d)
npz_file_basename = os.path.basename(npz_file_d)
npz_file_name = os.path.splitext(npz_file_basename)[0]
ID_1, ID_2 = npz_file_name.split("__vs__")

# Debugging (ctrl +1)
# print("Dirname:", npz_file_dirname)
# print("Basename:", npz_file_basename)
# print("Name:", npz_file_name)
# print("ID1:", ID_1)
# print("ID2:", ID_2)

def plot_coevolution(distance_matrix,
                     npz_file_name,
                     prefix,
                     a3m_seqs_N,
                     protein1="protein 1", protein2="protein 2",
                     title="RF2-track contacts"):

    # Plot the array as a heatmap
    plt.imshow(distance_matrix, cmap='hot', interpolation='nearest',
               
               # Max and minimum colors for heatmap
               vmin=0,
               vmax=1)
    
    # Labels
    plt.ylabel(protein1 + " (residue)")
    plt.xlabel(protein2 + " (residue)")
    plt.colorbar().set_label('Contact probability')
    
    # Add number of sequences
    plt.text(1.0, 1.10, f"Seqs Nº:\n{a3m_seqs_N}", ha='center',
                  va='center', transform=plt.gca().transAxes)
    
    # Add method that normalized the contact map (direct, reverse, min, max, mean)
    plt.text(0.0, 1.10, f"Method:\n{prefix}", ha='center',
              va='center', transform=plt.gca().transAxes)
        
    plt.title(title)

    # Save the plot to the specified output path
    output_path = prefix + "-" + npz_file_name + "-coevolution.png"
    plt.savefig(output_path)

    # Close the plot to release memory
    plt.close()

def obtain_metrics(distance_matrix,
                   top,
                   npz_file_name,
                   ID_1, ID_2,
                   prefix,
                   a3m_seqs_N):
    """
    

    Parameters
    ----------
    distance_matrix : array
        predicted contacts map array.
    top : int
        top contact residue pairs to extract.
    npz_file_dirname : str
        directory of the npz file.
    npz_file_name : str
        name of the npz file.
    ID_1 : str
        first protein ID.
    ID_2 : str
        second protein ID.
    prefix : str
        method that normalized the matrix (direct, reversed, minimum,
                                           maximum, mean).

    Returns
    -------
    metrics_df : pandas df
        contains several metrics for PPI prediction

    """
    
    # Flatten the array
    flat_matrix = distance_matrix.ravel()
    
    # Compute metrics to store as table
    metrics = {
        "ID1" : ID_1,
        "ID2" : ID_2,
        
        # Method (direct, reverse, mean, min, max)
        "method" : prefix,
                
        # Maximum contact probability value
        "max_contact" : distance_matrix.max(),
        
        # Sums greater than 0.0, 0.1, 0.2, etc
        "summ_000" : sum(flat_matrix),
        "summ_005" : sum(flat_matrix[flat_matrix > 0.05]),
        "summ_010" : sum(flat_matrix[flat_matrix > 0.10]),
        "summ_015" : sum(flat_matrix[flat_matrix > 0.15]),
        "summ_020" : sum(flat_matrix[flat_matrix > 0.20]),
        "summ_025" : sum(flat_matrix[flat_matrix > 0.25]),
        "summ_030" : sum(flat_matrix[flat_matrix > 0.30]),
        "summ_035" : sum(flat_matrix[flat_matrix > 0.35]),
        "summ_040" : sum(flat_matrix[flat_matrix > 0.40]),
        "summ_045" : sum(flat_matrix[flat_matrix > 0.45]),
        "summ_050" : sum(flat_matrix[flat_matrix > 0.50]),
        "summ_055" : sum(flat_matrix[flat_matrix > 0.55]),
        "summ_060" : sum(flat_matrix[flat_matrix > 0.60]),
        "summ_065" : sum(flat_matrix[flat_matrix > 0.65]),
        "summ_070" : sum(flat_matrix[flat_matrix > 0.70]),
        "summ_075" : sum(flat_matrix[flat_matrix > 0.75]),
        "summ_080" : sum(flat_matrix[flat_matrix > 0.80]),
        "summ_085" : sum(flat_matrix[flat_matrix > 0.85]),
        "summ_090" : sum(flat_matrix[flat_matrix > 0.90]),
        "summ_095" : sum(flat_matrix[flat_matrix > 0.95]),
        
        # Normalized sums
        "norm_summ_000" : sum(flat_matrix) / len(flat_matrix) if len(flat_matrix) != 0 else 0,
        "norm_summ_005" : sum(flat_matrix[flat_matrix > 0.05]) / len(flat_matrix[flat_matrix > 0.05]) if len(flat_matrix[flat_matrix > 0.05]) != 0 else 0,
        "norm_summ_010" : sum(flat_matrix[flat_matrix > 0.10]) / len(flat_matrix[flat_matrix > 0.10]) if len(flat_matrix[flat_matrix > 0.10]) != 0 else 0,
        "norm_summ_015" : sum(flat_matrix[flat_matrix > 0.15]) / len(flat_matrix[flat_matrix > 0.15]) if len(flat_matrix[flat_matrix > 0.15]) != 0 else 0,
        "norm_summ_020" : sum(flat_matrix[flat_matrix > 0.20]) / len(flat_matrix[flat_matrix > 0.20]) if len(flat_matrix[flat_matrix > 0.20]) != 0 else 0,
        "norm_summ_025" : sum(flat_matrix[flat_matrix > 0.25]) / len(flat_matrix[flat_matrix > 0.25]) if len(flat_matrix[flat_matrix > 0.25]) != 0 else 0,
        "norm_summ_030" : sum(flat_matrix[flat_matrix > 0.30]) / len(flat_matrix[flat_matrix > 0.30]) if len(flat_matrix[flat_matrix > 0.30]) != 0 else 0,
        "norm_summ_035" : sum(flat_matrix[flat_matrix > 0.35]) / len(flat_matrix[flat_matrix > 0.35]) if len(flat_matrix[flat_matrix > 0.35]) != 0 else 0,
        "norm_summ_040" : sum(flat_matrix[flat_matrix > 0.40]) / len(flat_matrix[flat_matrix > 0.40]) if len(flat_matrix[flat_matrix > 0.40]) != 0 else 0,
        "norm_summ_045" : sum(flat_matrix[flat_matrix > 0.45]) / len(flat_matrix[flat_matrix > 0.45]) if len(flat_matrix[flat_matrix > 0.45]) != 0 else 0,
        "norm_summ_050" : sum(flat_matrix[flat_matrix > 0.50]) / len(flat_matrix[flat_matrix > 0.50]) if len(flat_matrix[flat_matrix > 0.50]) != 0 else 0,
        "norm_summ_055" : sum(flat_matrix[flat_matrix > 0.55]) / len(flat_matrix[flat_matrix > 0.55]) if len(flat_matrix[flat_matrix > 0.55]) != 0 else 0,
        "norm_summ_060" : sum(flat_matrix[flat_matrix > 0.60]) / len(flat_matrix[flat_matrix > 0.60]) if len(flat_matrix[flat_matrix > 0.60]) != 0 else 0,
        "norm_summ_065" : sum(flat_matrix[flat_matrix > 0.65]) / len(flat_matrix[flat_matrix > 0.65]) if len(flat_matrix[flat_matrix > 0.65]) != 0 else 0,
        "norm_summ_070" : sum(flat_matrix[flat_matrix > 0.70]) / len(flat_matrix[flat_matrix > 0.70]) if len(flat_matrix[flat_matrix > 0.70]) != 0 else 0,
        "norm_summ_075" : sum(flat_matrix[flat_matrix > 0.75]) / len(flat_matrix[flat_matrix > 0.75]) if len(flat_matrix[flat_matrix > 0.75]) != 0 else 0,
        "norm_summ_080" : sum(flat_matrix[flat_matrix > 0.80]) / len(flat_matrix[flat_matrix > 0.80]) if len(flat_matrix[flat_matrix > 0.80]) != 0 else 0,
        "norm_summ_085" : sum(flat_matrix[flat_matrix > 0.85]) / len(flat_matrix[flat_matrix > 0.85]) if len(flat_matrix[flat_matrix > 0.85]) != 0 else 0,
        "norm_summ_090" : sum(flat_matrix[flat_matrix > 0.90]) / len(flat_matrix[flat_matrix > 0.90]) if len(flat_matrix[flat_matrix > 0.90]) != 0 else 0,
        "norm_summ_095" : sum(flat_matrix[flat_matrix > 0.95]) / len(flat_matrix[flat_matrix > 0.95]) if len(flat_matrix[flat_matrix > 0.95]) != 0 else 0,
        "paired_seqs_number" : a3m_seqs_N
        }
    
    # Create a pandas dataframe from the dictionary and save as tsv
    metrics_df = pd.DataFrame(metrics, index=[0])
    
    # Find the indices that would sort the "dist" array in descending order
    sorted_indices = np.argsort(distance_matrix, axis=None)[::-1]
    
    # Select the first <top> indices from the sorted list
    top_indices = sorted_indices[:top]
    
    # ChimeraX distances commands
    chimera_file = prefix + "-" + npz_file_name  + ".chimeraX"
    file_1 = open(chimera_file, "w")
    file_1.write("# ChimeraX commands to select the C-alpha of the predicted contacts" + '\n')
    
    # Top contacts file
    contacts_file = prefix + "-" + npz_file_name  + ".contacts"
    file_2 = open(contacts_file, "w")
    file_2.write(f"{ID_1}\t{ID_2}\tcontact_probability" + '\n')
    
    contacts = []
    for index in top_indices:
        
        # Store contacts as X and Y values
        x, y = np.unravel_index(index, distance_matrix.shape)
        
        # Contact probability value
        contact_prob = distance_matrix[x, y]
        
        # +1 to match the residue in visualization programs
        contacts.append([x+1, y+1, contact_prob])
        
        # Commands for ChimeraX
        file_1.write(f"distance /A:{x+1}@CA /B:{y+1}@CA" + '\n')
        file_2.write(f"{x+1}\t{y+1}\t{contact_prob}" + '\n')
    
    file_1.write("~label" + '\n')
    
    # Close the files
    file_1.close()
    file_2.close()
    
    return metrics_df


###############################################################################
######################### Preprocess contact maps #############################
###############################################################################

# Load the npz files and extract the "dist" array -----------------------------

# Direct data
data_d = np.load(npz_file_d)
dist_d = data_d['dist']

# Reversed transposed data
data_r = np.load(npz_file_r)
dist_r = data_r['dist']
dist_rt =  np.transpose(dist_r)

# Generate derivates ----------------------------------------------------------

# Minumum
dist_min = np.minimum(dist_d, dist_rt)

# Maximum
dist_max = np.maximum(dist_d, dist_rt)

# Mean
dist_mean = (dist_d + dist_rt) / 2


###############################################################################
######################### Find predicted contacts #############################
###############################################################################

# Obtain the metrics of each normalization method
mertics_df_d = obtain_metrics(dist_d, top, npz_file_name, ID_1, ID_2, "direct", a3m_seqs_N)
mertics_df_rt = obtain_metrics(dist_rt, top, npz_file_name, ID_1, ID_2, "reversed", a3m_seqs_N)
mertics_df_min = obtain_metrics(dist_min, top, npz_file_name, ID_1, ID_2, "min", a3m_seqs_N)
mertics_df_max = obtain_metrics(dist_max, top, npz_file_name, ID_1, ID_2, "max", a3m_seqs_N)
mertics_df_mean = obtain_metrics(dist_mean, top, npz_file_name, ID_1, ID_2, "mean", a3m_seqs_N)

# Merge the metrics dataframes
full_metrics_df = pd.concat([mertics_df_d,
                             mertics_df_rt,
                             mertics_df_min,
                             mertics_df_max,
                             mertics_df_mean],
                            # Merge rows axis
                            axis=0)

# Save DF as tsv
metrics_file = npz_file_name + ".metrics"
full_metrics_df.to_csv(metrics_file , sep='\t', index=False)

###############################################################################
######################## Plot coevolution heatmap #############################
###############################################################################

plot_coevolution(dist_d,
                 npz_file_name,
                 "direct",
                 a3m_seqs_N,
                 protein1 = ID_1,
                 protein2 = ID_2,
                 title="RF2-track contacts")
plot_coevolution(dist_rt,
                 npz_file_name,
                 "reversed",
                 a3m_seqs_N,
                 protein1 = ID_1,
                 protein2 = ID_2,
                 title="RF2-track contacts")
plot_coevolution(dist_min,
                 npz_file_name,
                 "min",
                 a3m_seqs_N,
                 protein1 = ID_1,
                 protein2 = ID_2,
                 title="RF2-track contacts")
plot_coevolution(dist_max,
                 npz_file_name,
                 "max",
                 a3m_seqs_N,
                 protein1 = ID_1,
                 protein2 = ID_2,
                 title="RF2-track contacts")
plot_coevolution(dist_mean,
                 npz_file_name,
                 "mean",
                 a3m_seqs_N,
                 protein1 = ID_1,
                 protein2 = ID_2,
                 title="RF2-track contacts")
