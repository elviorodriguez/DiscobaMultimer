#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 12:18:04 2023

@author: elvio

This script performs the pairing of two a3m MSA
"""

import sys
import os
from Bio import SeqIO
import re
from Bio import pairwise2
from Bio.Align import PairwiseAligner

###############################################################################
############################## Checking and usage #############################
###############################################################################

# Check that two command-line arguments have been provided
if len(sys.argv) < 4:
    print("Error: At least two a3m files and an output name are required", file=sys.stderr)
    print("USAGE: python perform_pairing.py <output.a3m> <prot_ID1.a3m> <prot_ID2.a3m> [<prot_IDn.a3m>]", file=sys.stderr)
    print("   output.a3m        name and path of the outputh a3m file to save")
    print("   prot_ID1.a3m      first .a3m file generated with run_MMseqs2_to_get_DiscobaMSA.sh")
    print("   prot_ID2.a3m      second .a3m file generated with run_MMseqs2_to_get_DiscobaMSA.sh")
    print("   prot_IDn.a3m      Nth .a3m file (optional and as many as needed)")
    print("OUTPUT:")
    print("   output.a3m     : paired+unpaired a3m file with cardinality")
    sys.exit(1)

# Assign the command-line arguments to variables
output_file = sys.argv[1]   # output file name and path
a3m_files={}
for i in range(2, len(sys.argv)):
    a3m_files[f"a3m_{i-1}"] = sys.argv[i]

# # Debug
# print(output_file)
# print(a3m_files.keys())
# print(a3m_files.values())

# Check if input files are in .a3m format
for a3m_file in a3m_files.values():
    if not a3m_file.endswith('.a3m'):
        print("ERROR: Input files must be in .a3m format.", file=sys.stderr)
        sys.exit(1)

###############################################################################
############################### Helper functions ##############################
###############################################################################

# Function to check if it is an homooligomer
def is_homooligomer(a3m_files):
    """
    Checks if the prediction is comming from an homooligomer

    Parameters
    ----------
    a3m_files : dict
        Defined at the begining.

    Returns
    -------
    Bool.

    """
    subunits = len(a3m_files.keys())
    if len(set(a3m_files.values())) == 1:
        print(f"HOMO-oligomer with {subunits} subunits")
        return True
    else:
        print(f"HETERO-oligomer with {subunits} subunits")
        return False

def pair_homooligomer(a3m_files, output_file):
    
    subunits = len(a3m_files.keys())
    file_location=a3m_files["a3m_1"]
    with open(file_location, "r") as file_read:
       L = len(file_read.readlines()[1].rstrip("\n"))   # protein length
    
    # Cardinality line (e.g. #192,125	1,1)
    Ls=','.join(str(i) for i in [L]*subunits)
    Ns=','.join(['1']*subunits)
    cardinality = f"#{Ls}\t{Ns}"
    
    # Open the file in read mode
    with open(file_location, "r") as file_read:
        
        with open(output_file, "w") as file_write:
            
            # Add cardinality line
            file_write.write(cardinality + "\n")
            
            # Pair the sequences
            for line in file_read:
                if line.startswith(">"):
                    paired_header = line.rstrip("\n").lstrip(">")
                    file_write.write(">" + '\t'.join([paired_header] * subunits) + "\n")
                else:
                    paired_seq = line.rstrip("\n") * subunits
                    file_write.write(paired_seq + "\n")
                    

    
# Function to parse the a3m_file and group the sequences by TaxID in a dictionary
def separate_by_tax_id(a3m_file):
    """
    Parse the a3m_file and group the sequences by TaxID in a dictionary.
    
    Parameters
    ----------
    a3m_file : a3m file
        An a3m file generated with run_MMseqs2_to_get_DiscobaMSA.sh

    Returns
    -------
    records_by_taxid : dict
        Discoba TaxIDs as keys.
        List of sequence records as values, grouped by TaxID

    """
    
    # Load the a3m files
    records = SeqIO.parse(a3m_file, "fasta")
    
    # Create an empty dictionary to store the records in the first a3m file
    records_by_taxid = {}
    
    # Parse 
    for i, record in enumerate(records):
        
        # If it is the first record
        if i == 0:
            # Add "query" as TaxID annotation to the first sequence
            tax_id="query"
            record.annotations["TaxID"] = tax_id
            
        # For the rest, extract the TaxID and add it as annotation
        else:
            # Search for the TaxID pattern in the ID string
            match = re.search(r"TaxID=(\w+)", record.description)
            
            # # Debug
            # print(record.description)
            # print("match:", match)
            # print("")
            
            # print("tax_id:", match.group(1))
            if match:
                # Extract the TaxID from the first matching group
                tax_id = match.group(1)
                # Add the TaxID as an annotation
                record.annotations["TaxID"] = tax_id
        
        # Group the records by TaxID in a dictionary
        if tax_id not in records_by_taxid:  # Add new key if it does not exist
            records_by_taxid[tax_id] = []
        records_by_taxid[tax_id].append(record)
        
    # Return the dictionary of records grouped by TaxID
    return records_by_taxid

# # Debug
# print(separate_by_tax_id(a3m_files["a3m_1"]).keys())
# print("OK")

# Define a function to perform global pairwise alignment
def global_alignment(query_seq, subject_seq):
    """
    Perform a global pairwise alignment between two sequences using the Needleman-Wunsch algorithm.
    
    Parameters
    ----------
    query_seq : str
        The query sequence.
    subject_seq : str
        The subject sequence to be aligned to the query.
        
    Returns
    -------
    Tuple containing:
        - Aligned query sequence
        - Aligned subject sequence
        - Alignment score
    """
    
    alignments = pairwise2.align.globalxx(query_seq, subject_seq)
    best_alignment = alignments[0]  # Get the highest-scoring alignment
    aligned_query = best_alignment[0]
    aligned_subject = best_alignment[1].replace('-', '')
    alignment_score = best_alignment[2]
    
    coverage = len(aligned_subject) / len(query_seq) * 100
    
    return aligned_query, aligned_subject, alignment_score, coverage


def get_N_seq_by_TaxID(taxid_grouped):
    for TaxID in taxid_grouped:
        print(TaxID, "=", len(taxid_grouped[TaxID]))

def get_all_keys(dicts):
    """
    Takes a list of dictionaries and 

    Parameters
    ----------
    *dicts : list of dicts
        Any dictionary.

    Returns
    -------
    list
        Contains all the possible key in the dictionaries (without repetition).

    """
    keys = set()
    for d in dicts:
        keys.update(d.keys())
    return list(keys)        

def sort_by_similarity_to_query(taxid_grouped):
    """
    Sorts the sequence records inside the TaxID groups from highest to lowest
    similarity to the query and modifies the input dictionay.
    Adds the annotation "similarity_to_query"
    
    Parameters
    ----------
    taxid_grouped : dict
        DESCRIPTION.

    Returns
    -------
    Nothing
    
        It modify the input dict by sorting the sequences in each TaxID from
        highest to lowest similarity to the query and adds the annotation
        "similarity_to_query" to each sequence record (except for query).

    """
   
    # Iterate over each TaxID group and compute the similarity to the query
    for TaxID in taxid_grouped.keys():
        
        # Skip calculations for the query sequence
        if TaxID == "query":
            taxid_grouped["query"][0].annotations["similarity_to_query"] = 0
            continue

        # Iterate over each seq in the TaxID group
        for record in taxid_grouped[TaxID]:
            
            # Perform global alignment with the query and extract the score
            aligned_query, aligned_subject, alignment_score, coverage = global_alignment(
                query_seq = taxid_grouped["query"][0].seq,
                subject_seq = record.seq)
            
            record.annotations["similarity_to_query"] = alignment_score
            record.annotations["coverage_to_query"] = coverage
                
            # Sort the records in the TaxID group by their similarity to the query
            taxid_grouped[TaxID] = sorted(taxid_grouped[TaxID],
                                          key=lambda x: x.annotations.get("similarity_to_query", 0),
                                          # From highest to lowest
                                          reverse=True)

# # Debug
# grouped_by_taxid_a3m = separate_by_tax_id(a3m_files["a3m_1"])
# print(grouped_by_taxid_a3m)
# sort_by_similarity_to_query(grouped_by_taxid_a3m)
# print(grouped_by_taxid_a3m)
# print("OK")
# exit()             

def sort_protein_sequences_by_similarity(input_file, output_file):
    # read in fasta file
    records = list(SeqIO.parse(input_file, "fasta"))
    # set first sequence as reference
    reference_seq = records[0].seq
    # create PairwiseAligner object for global alignment
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    # calculate similarity scores for each sequence
    similarity_scores = []
    for record in records[1:]:
        alignment = aligner.align(reference_seq, record.seq)
        similarity_score = alignment.score / max(len(reference_seq), len(record.seq))
        similarity_scores.append((record, similarity_score))
    # sort sequences by decreasing similarity
    sorted_seqs = [records[0]]
    for seq, score in sorted(similarity_scores, key=lambda x: x[1], reverse=True):
        sorted_seqs.append(seq)
    # write sorted sequences to output file without newlines
    with open(output_file, "w") as handle:
        for seq in sorted_seqs:
            handle.write(">" + seq.description + "\n" + str(seq.seq) + "\n")
    return

def insert_string_as_first_line(file_path, string):
    # read in contents of file
    with open(file_path, "r") as f:
        contents = f.read()
    # prepend string to contents
    updated_contents = string + "\n" + contents
    # write updated contents back to file
    with open(file_path, "w") as f:
        f.write(updated_contents)



def pair_a3m(grouped_sorted_a3m_list, output_file, min_coverage):
    """
    Performs pairing alignment and generates a paired.a3m file.

    Parameters
    ----------
    grouped_sorted_a3m_files : list (with dicts)
        Contains at least two dict with the grouped by TaxID sequences sorted
        by similarity to the query.
        
    output_file : str
        Expected to be in this format: "<protID_1>__vs__<protID_2>.a3m".
        E.g.: "C4B63_40g26__vs__C4B63_28g81.a3m"
        If more that two query proteins are used, they will be all separated by
        "__vs__".
        E.g.:"C4B63_40g26__vs__C4B63_40g26__vs__C4B63_28g81.a3m"

    Returns
    -------
    None. Just generates the output_file

    """
    
    # Remove query sequence records
    queries = {"records" : [],
               "IDs" : [],
               "Lengths" : [],
               "Seqs" : []}
    
    for i, grouped_sorted_a3m in enumerate(grouped_sorted_a3m_list):
        
        # Query record
        record = grouped_sorted_a3m.pop("query")[0]
        queries["records"].append(record)
    
        # Query IDs
        queries["IDs"].append(record.id)
    
        # Query Lengths
        queries["Lengths"].append(len(record))
    
        # Query Sequences
        queries["Seqs"].append(str(record.seq))
        
    
    # Cardinality line (e.g. #192,125	1,1)
    Ls=','.join(str(i) for i in queries["Lengths"])
    Ns=','.join(['1']*len(queries["Lengths"]))
    IDs='\t'.join(queries["IDs"])
    cardinality = f"#{Ls}\t{Ns}"
    paired_query_header=f">{IDs}"
    paired_query_seq=''.join(queries["Seqs"])
 
    # # Debug
    # print(queries)
    # print(all_TaxIDs)
    # print(cardinality)
    # print(paired_query_header)
    # print(paired_query_seq)
            
    # List all TaxIDs covered by all proteins
    all_TaxIDs = get_all_keys(grouped_sorted_a3m_list)
        
    # Open a file for writing
    with open(output_file, 'w') as f:
        
        # Write cardinality and paired queries      
        f.write(cardinality+'\n')
        f.write(paired_query_header+'\n')
        f.write(paired_query_seq+'\n')
        
        # Parse the sequences by groups
        for TaxID in all_TaxIDs:
            
            # Check if all have sequences for the TaxID
            if not all(TaxID in d for d in grouped_sorted_a3m_list):
                # Debug
                # print(f'The key "{TaxID}" is not in one of the dictionaries')
                continue
            
            else:
                # Debug
                # print(f'TaxID: {TaxID} ----------------------------------------')
                
                # See which query has more sequences for the same TaxID (notice that if any of the query sequences did not retrived an homolog at an specific TaxID, the paired alignment for that TaxID will be omitted)
                min_N_seq = min(len(position[TaxID]) for position in grouped_sorted_a3m_list)
                
                # The list of seq for the TaxID is sorted by similarity to query
                for rank in range(min_N_seq):
                    
                    # Check if all have at least min_coverage % to the query
                    all_covered = True
                    for position in grouped_sorted_a3m_list:
                        if position[TaxID][rank].annotations["coverage_to_query"] < min_coverage:
                            all_covered = False
                            break
                    
                    if all_covered:
                        # Store subject IDs and sequences
                        subjects = {"IDs" : [],
                                   "Seqs" : []}
                                                               
                        # For each a3m position
                        for position in grouped_sorted_a3m_list:
                            # Subject IDs
                            subjects["IDs"].append(position[TaxID][rank].id)
                            
                            # Subject Seqs
                            subjects["Seqs"].append(str(position[TaxID][rank].seq))
                            
                        # Pair the result
                        paired_subj_header=">" + '\t'.join(subjects["IDs"])
                        paired_subj_seq=''.join(subjects["Seqs"])
                        
                        # # debug
                        # print(paired_subj_header)
                        # print(paired_subj_seq)
                            
                        # Write paired
                        f.write(paired_subj_header+'\n')
                        f.write(paired_subj_seq+'\n')
    
    sort_protein_sequences_by_similarity(output_file, "temporal_file")
    insert_string_as_first_line("temporal_file", cardinality)
    os.rename("temporal_file", output_file)
    
# # Debug
# grouped_by_taxid_a3m = separate_by_tax_id(a3m_files["a3m_1"])
# sort_by_similarity_to_query(grouped_by_taxid_a3m)
# copy = grouped_by_taxid_a3m.copy()
# copy2 = grouped_by_taxid_a3m.copy()
# grouped_sorted_a3m_list = [grouped_by_taxid_a3m, copy, copy2]

# pair_a3m(grouped_sorted_a3m_list, output_file)
# print("OK")
# exit()

def add_unpaired(output_file, a3m_files):
    """
    Appends the unpaired part to the paired a3m file generated by pair_a3m().

    Parameters
    ----------
    output_file : str
        Path to the file generated by pair_a3m()
        Expected to be in this format "<protID_1>__vs__<protID_2>.a3m".
        E.g.: "C4B63_40g26__vs__C4B63_28g81.a3m".
    a3m_files : dict
        Contains the paths to the a3m files for the succesive paired sequence.
        
    Returns
    -------
    None. It only appends the unpaired alignment to the output_file.

    """
    
    queries = {"IDs" : [],
               "Seqs" : [],
               "Lengths" : []}
    
    # Extract query sequence info
    for a3m_position in a3m_files.keys():
        
        # Load records
        records = SeqIO.parse(a3m_files[a3m_position], "fasta")
        
        # Extract queries
        query = next(records)
    
        # Queries seqs
        query_header = query.id
        
        # Queries seqs
        query_seq = str(query.seq)
        
        # Queries lengths
        L = len(query.seq)
        
        # Add data to dict
        queries["IDs"].append(query_header)
        queries["Seqs"].append(query_seq)
        queries["Lengths"].append(L)
    
    # Number of queries
    N_positions = len(queries["IDs"])
        
    # Make unpaired part
    for a3m_key in a3m_files.keys():
        
        # Cero based position
        position=int(a3m_key.split('_')[1])-1
        
        # Load records
        records = SeqIO.parse(a3m_files[a3m_key], "fasta")
        
        # Jump query
        next(records)
        
        # Open a file for writing
        with open(output_file, 'a') as f:
            
            # Fill with gaps
            query_seq = []
            for i in range(N_positions):
                if i == position:
                    query_seq.append(queries["Seqs"][i])
                else:
                    query_seq.append("-" * queries["Lengths"][i])
                        
            # Write header for the query
            f.write(">" + queries["IDs"][position] + '\n')
            f.write(''.join(query_seq) + '\n')
        
            # Add unpaired subjects frome the first a3m file
            for record in records:
                
                # Extract ID and sequence
                subj_ID = record.id
                subj_seq = str(record.seq)
                
                # Fill with gaps
                gap_subj_seq = []
                for i in range(N_positions):
                    if i == position:
                        gap_subj_seq.append(subj_seq)
                    else:
                        gap_subj_seq.append("-" * queries["Lengths"][i])
                
                # Write them with pad gaps
                f.write(">" + subj_ID + '\n')
                f.write(''.join(gap_subj_seq) + '\n')

        
# # Debug
# grouped_by_taxid_a3m = separate_by_tax_id(a3m_files["a3m_1"])
# sort_by_similarity_to_query(grouped_by_taxid_a3m)
# copy = grouped_by_taxid_a3m.copy()
# copy2 = grouped_by_taxid_a3m.copy()
# grouped_sorted_a3m_list = [grouped_by_taxid_a3m, copy, copy2]

# pair_a3m(grouped_sorted_a3m_list, output_file)

# add_unpaired(output_file, a3m_files)
# print("OK")
# exit()

###############################################################################
################################### Running ###################################
###############################################################################

# Set min coverage for pairing at 50%
min_coverage = 50 

# Group a3m files sequences by TaxID and rank them by similarity to query
a3m_taxid =  {}
for a3m in a3m_files.keys():
    a3m_taxid[a3m] = separate_by_tax_id(a3m_files[a3m])
    sort_by_similarity_to_query(a3m_taxid[a3m])
    
# Convert a3m_taxid to list of dict
grouped_sorted_a3m_list = []
for a3m in a3m_taxid.keys():
    grouped_sorted_a3m_list.append(a3m_taxid[a3m])
    
# Generate a file with the paired part
pair_a3m(grouped_sorted_a3m_list, output_file, min_coverage)
    
# Append the unpaired part to the file
add_unpaired(output_file, a3m_files)

# Check if sorting by similarity to query has worked
# for TaxID in taxid_1.keys():
#     # Skip the query sequence
#     if TaxID == "query":
#         continue
#     print("TaxID:", TaxID, "--------------------------------------------------")
#     print()
#     for i in range(len(taxid_1[TaxID])):
#         print("     ", taxid_1[TaxID][i].id)
#         # print("         Seq:", taxid_1[TaxID][i].seq)
#         print("         Rank:", i)
#         print("         Similarity to query:", taxid_1[TaxID][i].annotations["similarity_to_query"])
#         print()



