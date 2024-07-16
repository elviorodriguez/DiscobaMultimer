#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 12:18:04 2023

@author: elvio

This script performs the pairing of two or more a3m discoba MSAs
"""
# Remove deprecation warning
import warnings
from Bio import BiopythonDeprecationWarning
warnings.simplefilter(action='ignore', category=BiopythonDeprecationWarning)

import sys
import os
from Bio import SeqIO, pairwise2
import re
from Bio.Align import PairwiseAligner
from collections import Counter

###############################################################################
############################## Checking and usage #############################
###############################################################################

# Check that two command-line arguments have been provided
if len(sys.argv) < 4:
    print("Error: At least two a3m files and an output name are required", file=sys.stderr)
    print("USAGE: python perform_greedy_pairing.py <output.a3m> <prot_ID1.a3m> <prot_ID2.a3m> [<prot_IDn.a3m>]", file=sys.stderr)
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
# print("split:", str(perform_splitting))

# Check if input files are in .a3m format
for a3m_file in a3m_files.values():
    if not a3m_file.endswith('.a3m'):
        print("ERROR: Input files must be in .a3m format.", file=sys.stderr)
        sys.exit(1)

###############################################################################
############################### Helper functions ##############################
###############################################################################

def index_exists(lst, index):
    try:
        lst[index]
        return True
    except IndexError:
        return False

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


def is_protein_repeated(a3m_files):
    """
    Returns true if the prediction has a protein repeated more than once.

    Parameters
    ----------
    a3m_files : dict
        Defined at the begining.

    Returns
    -------
    Bool.

    """
    subunits = len(a3m_files.values())
    subunits_set = len(set(a3m_files.values()))
    
    # Compare the length of the set with the length of the original list
    if subunits_set < subunits:
        return True  # Repeated proteins found
    else:
        return False  # No repeated proteins

# print(is_protein_repeated(a3m_files))


def find_number_of_each_proteins(a3m_files):
    """
    Returns the number of each protein 

    Parameters
    ----------
    a3m_files : dic
        Defined at the begining.

    Returns
    -------
    element_counts : list
        A list of tuples with the number of times each protein (a3m file)
        appears.
        eg: [('Tb10.v4.0245.a3m', 3), ('Tb927.3.5740.a3m', 1)]

    """
    counter = Counter(a3m_files.values())
    proteins_counts = counter.most_common()
    return proteins_counts

# IMPLEMENTED IN discoba-multimer_batch.sh (not used here)
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
                    
# IMPLEMENTED IN discoba-multimer_batch.sh (not used here)
def pair_homooligomer_NOsplit(a3m_files, output_file):
    """
    Creates an a3m file (output_file) with a single sequence unpaired, adding
    the cardinality with the number of subunits in the homooligomer
    eg: "#102	2" (protein length 102 and 2 subunits)

    Parameters
    ----------
    a3m_files : dic
        Defined at the begining.

    output_file : str
        output file name (.a3m).

    Returns
    -------
    None.

    """
    
    subunits = len(a3m_files.keys())
    file_location=a3m_files["a3m_1"]
    
    with open(file_location, "r") as file_read:
       L = len(file_read.readlines()[1].rstrip("\n"))   # protein length
    
    # Cardinality line (e.g. #192	2)
    Ls=L
    Ns=subunits
    cardinality = f"#{Ls}\t{Ns}"
    
    # Open the file in read mode
    with open(file_location, "r") as file_read:
        
        with open(output_file, "w") as file_write:
            
            # Add cardinality line
            file_write.write(cardinality + "\n")
            
            # Add the sequences
            for line in file_read:
                file_write.write(line)

# pair_homooligomer_NOsplit(a3m_files, "test.a3m")
# sys.exit()
    
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
    
    # Remove gaps
    query_seq = query_seq.replace('-', '')
    subject_seq = subject_seq.replace('-', '')
    
    # Compute 
    alignments = pairwise2.align.globalxx(query_seq, subject_seq)
    best_alignment = alignments[0]  # Get the highest-scoring alignment
    aligned_query = best_alignment.seqA
    aligned_subject = best_alignment.seqB
    alignment_score = best_alignment.score
    
    # Check if something went wrong
    if len(aligned_query) != len(aligned_subject): 
        raise ValueError(f"FATAL ERROR: something went wrong while performing global alignments of the following sequences: \nSEQUENCE_1: {query_seq}\nSEQUENCE_2: {subject_seq}")
    
    # Compute sequence coverage
    match_counter = 0
    for pos in range(len(aligned_query)):
        if aligned_query[pos] != "-" and aligned_subject[pos] != "-":
            match_counter += 1
    coverage = match_counter / len(query_seq) * 100
    
    return aligned_query, aligned_subject, alignment_score, coverage

# # Debugging
# query_seq   = "METTLASYPLIRASTREWQYTVCMNH".replace('-', '')
# subject_seq = "qqMEaaaTTLASYPLIRaaaASTREWQYaaaTVCMNH".replace('-', '')

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

def add_similarity_to_query(taxid_grouped):
    """
    Adds the annotations "similarity_to_query" and "coverage_to_query" to each
    sequence record.
    
    Parameters
    ----------
    taxid_grouped : dict
        Dictionary containing TaxID as keys and sequences as values.

    Returns
    -------
    Nothing
    
        It modify the input dict by sorting the sequences in each TaxID from
        highest to lowest similarity to the query and adds the annotation
        "similarity_to_query" and "coverage_to_query" to each sequence record 
        as annotation.

    """
   
    # Iterate over each TaxID group and compute the similarity to the query
    for TaxID in taxid_grouped.keys():
        
        # Skip calculations for the query sequence
        if TaxID == "query":
            taxid_grouped["query"][0].annotations["similarity_to_query"] = 0
            taxid_grouped["query"][0].annotations["coverage_to_query"] = 100
            continue

        # Iterate over each seq in the TaxID group
        for record in taxid_grouped[TaxID]:
            
            # Perform global alignment with the query and extract the score
            aligned_query, aligned_subject, alignment_score, coverage = global_alignment(
                query_seq = str(taxid_grouped["query"][0].seq),
                subject_seq = str(record.seq))
            
            record.annotations["similarity_to_query"] = alignment_score
            record.annotations["coverage_to_query"] = coverage
                
def remove_low_coverage_sequences(taxid_grouped_annotated, min_coverage):
    """
    Removes the low coverage sequences (< min_coverage) from TaxID dict

    Parameters
    ----------
    taxid_grouped_annotated : dict
        Dictionary containing TaxID as keys and sequences as values. Each
        sequence record contains "similarity_to_query" and "coverage_to_query"
        as annotations.
    min_coverage : flot (0 to 100)
        Cutoff coverage value

    Returns
    -------
    None. Modifies the dictionaries by removing sequences with low coverage

    """
    for TaxID in list(taxid_grouped_annotated.keys()):
        filtered_sequences = [
            seq for seq in taxid_grouped_annotated[TaxID]
            if seq.annotations["coverage_to_query"] >= min_coverage
        ]
        
        if len(filtered_sequences) != len(taxid_grouped_annotated[TaxID]):
            for seq in taxid_grouped_annotated[TaxID]:
                if seq.annotations["coverage_to_query"] < min_coverage:
                    print(f"Removing sequence: {seq}")
                    print(f"Coverage: {seq.annotations['coverage_to_query']} < {min_coverage}")

        taxid_grouped_annotated[TaxID] = filtered_sequences
 

def sort_by_similarity_to_query(taxid_grouped_annotated):
    """
    Sorts the sequence records inside the TaxID groups from highest to lowest
    similarity to the query and modifies the input dictionay.
    
    Parameters
    ----------
    taxid_grouped : dict
        DESCRIPTION.

    Returns
    -------
    Nothing
    
        It modify the input dict by sorting the sequences in each TaxID from
        highest to lowest similarity to the query (except for query).

    """
   
    # Iterate over each TaxID group and compute the similarity to the query
    for TaxID in taxid_grouped_annotated.keys():
        
        # Skip calculations for the query sequence
        if TaxID == "query":
            continue

        # Iterate over each seq in the TaxID group
        for record in taxid_grouped_annotated[TaxID]:

            # Sort the records in the TaxID group by their similarity to the query
            taxid_grouped_annotated[TaxID] = sorted(
                taxid_grouped_annotated[TaxID],
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
        alignment = aligner.align(str(reference_seq).replace("-", ""),
                                      str(record.seq).replace("-", ""))

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



def pair_a3m(grouped_sorted_a3m_list, output_file):
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

    # Paired queries as fasta (first paired sequeces)
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
        
    # Open output file for writing
    with open(output_file, 'w') as f:
        
        # Write cardinality and paired queries (first 3 lines)
        f.write(cardinality+'\n')
        f.write(paired_query_header+'\n')
        f.write(paired_query_seq+'\n')
        
        # Parse the sequences by TaxID group
        for TaxID in all_TaxIDs:
            
            # Skip pairing if TaxID is in at least two sequences
            if sum(TaxID in d for d in grouped_sorted_a3m_list) < 2:
                # # Debug
                # print(f'The key "{TaxID}" is not in at least two of the dictionaries')
                continue

            else:
                # # Debug
                # print(f'TaxID: {TaxID} ----------------------------------------')
                                    
                # See which is the query with more sequences for the same TaxID
                max_N_seq = max(len(position.get(TaxID, [])) for position in grouped_sorted_a3m_list)
                
                # Make the pairing one rank at a time (Seqs are sorted by similarity to query)
                for rank in range(max_N_seq):
                    
                    # Boolean list indicating which sequence have TaxID at current rank (eg [True, False, True])
                    bool_list = list(index_exists(position.get(TaxID, []), rank) for position in grouped_sorted_a3m_list)
                    
                    # Skip pairing if TaxID is not in at least two sequences for the rank
                    if sum(bool_list) < 2:
                        # # Debug
                        # print(f'The key "{TaxID}" is not in at least two of the dictionaries for the rank {rank}')
                        continue
                    
                    # Store subject IDs and sequences
                    subjects = {"IDs" : [],
                               "Seqs" : []}
                                                           
                    # For each a3m position
                    for i, position in enumerate(grouped_sorted_a3m_list):
                        
                        # If the position have a protein in the rank
                        if bool_list[i]:
                            # Subject IDs
                            subjects["IDs"].append(position[TaxID][rank].id)
                            
                            # Subject Seqs
                            subjects["Seqs"].append(str(position[TaxID][rank].seq))
                        
                        # If the position does not have a protein at a given rank
                        else:
                            # Add empty ID
                            subjects["IDs"].append(f"no_rank_{rank}")
                            
                            # Fill with gaps
                            subjects["Seqs"].append(str(int(queries["Lengths"][i])*"-"))
                        
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

# Removes duplicated a3m files from the dictionary
def remove_duplicates(dictionary):
    output_dict = {}
    values_seen = set()

    for key, value in dictionary.items():
        if value not in values_seen:
            output_dict[key] = value
            values_seen.add(value)

    return output_dict

def update_dictionary_keys(dictionary):
    new_dict = {}
    key_count = 1

    for key, value in dictionary.items():
        new_key = 'a3m_' + str(key_count)
        new_dict[new_key] = value
        key_count += 1

    return new_dict


###############################################################################
################################### Running ###################################
###############################################################################

# # DEBUGGING input -------------------------------------------------------------

# def create_concatenated_filename(file_paths):
#     """
#     Takes a list of file paths and returns a string where the base names of the files
#     are concatenated with '__vs__' and the '.a3m' extension is added at the end.

#     :param file_paths: List of file paths
#     :return: Concatenated string
#     """
#     # Extract the base names without extensions
#     base_names = [os.path.splitext(os.path.basename(path))[0] for path in file_paths]
    
#     # Join the names with '__vs__' separator
#     joined_names = '__vs__'.join(base_names)
    
#     # Add the .a3m extension
#     result = f"{joined_names}.a3m"
    
#     return result


# os.chdir('/home/elvio/DiscobaMultimer/scripts')
# input_a3ms_1 = ["../testing/BDF2.a3m", "../testing/TIF2.a3m"]
# input_a3ms_2 = ["../testing/BDF2.a3m", "../testing/TelAP1.a3m", "../testing/TelAP2.a3m"]
# input_a3ms_3 = ["../testing/BDF2.a3m", "../testing/TIF2.a3m", "../testing/TRF.a3m"]

# to_test = input_a3ms_1

# output_file = "../testing/" + create_concatenated_filename(to_test)

# a3m_files={}
# for i, a3m_f in enumerate(to_test):
#     a3m_files[f"a3m_{i+1}"] = a3m_f

# # # Debug
# # print(output_file)
# # print(a3m_files.keys())
# # print(a3m_files.values())

# # Check if input files are in .a3m format
# for a3m_file in a3m_files.values():
#     if not a3m_file.endswith('.a3m'):
#         print("ERROR: Input files must be in .a3m format.", file=sys.stderr)
#         sys.exit(1)


# Running ---------------------------------------------------------------------


# Set min coverage for pairing at 10%
min_coverage = 10

# Check if there is any protein repeated in the list
is_repeated = is_protein_repeated(a3m_files)
# If there is any repeated protein, remove duplicates and keep track of subunit
# number
if is_repeated:
    each_protein_number = find_number_of_each_proteins(a3m_files)
    each_protein_number = sorted(each_protein_number, key=lambda x: list(a3m_files.values()).index(x[0]))
    a3m_files = update_dictionary_keys(remove_duplicates(a3m_files))

# Group a3m files sequences by TaxID, remove low coverage and sort
a3m_taxid =  {}
for a3m in a3m_files.keys():
    
    # Grouping by TaxID
    a3m_taxid[a3m] = separate_by_tax_id(a3m_files[a3m])
    
    # Annotate
    add_similarity_to_query(taxid_grouped = a3m_taxid[a3m])
    

    # Remove sequences with less than min_coverage % to query
    remove_low_coverage_sequences(taxid_grouped_annotated = a3m_taxid[a3m],
                                  min_coverage = min_coverage)
    
    # Rank them by similarity to query
    sort_by_similarity_to_query(a3m_taxid[a3m])


# Convert a3m_taxid to list of dict
grouped_sorted_a3m_list = []
for a3m in a3m_taxid.keys():
    grouped_sorted_a3m_list.append(a3m_taxid[a3m])

# Generate a file with the paired part
pair_a3m(grouped_sorted_a3m_list, output_file)
    
# Append the unpaired part to the file
add_unpaired(output_file, a3m_files)

# Renumber the cardinality to match the subunit numbers
if is_repeated:
    
    # Read the contents of the file
    with open(output_file, 'r') as file:
        lines = file.readlines()
    
    # Modify the first line
    old_cardinality = lines[0].split('\t')
    new_cardinality = old_cardinality[0] + '\t' + ','.join([str(i[1]) for i in each_protein_number]) + '\n'
    lines[0] = new_cardinality
    
    # # Debug
    # print("old:", old_cardinality)
    # print("new:", new_cardinality)

    # Write the modified contents back to the file
    with open(output_file, 'w') as file:
        file.writelines(lines)


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