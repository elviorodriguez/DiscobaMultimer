#!/usr/bin/env python3

import itertools
import sys
import os

def generate_combinations(input_file, output_file, exclude_self):
    try:
        # Read input file
        with open(input_file, 'r') as file:
            ids = [line.strip() for line in file.readlines()]

        # Generate combinations
        combinations = list(itertools.combinations_with_replacement(ids, 2))
        
        # If -e flag was passed
        if exclude_self:
            combinations = [comb for comb in combinations if comb[0] != comb[1]]

        # Write output file
        with open(output_file, 'w') as file:
            for combo in combinations:
                file.write('\t'.join(combo) + '\n')

        print(f"Combinations generated successfully and saved to {output_file}")

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

def usage():
    print("Usage: ./generate_2mers_combinations.py <single_IDs.txt> <output_IDs_table.txt> [-e]")
    print("")
    print("Parameters:")
    print("   single_IDs.txt        : File with one ID/protein name per line.")
    print("   output_IDs_table.txt  : Name of the output TSV file.")
    print("   -e (optional)         : Excludes homodimers combinations.")
    sys.exit(1)

if __name__ == "__main__":
    # Check if correct number of arguments is provided
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        usage()

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    exclude_self = "-e" in sys.argv

    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)

    generate_combinations(input_file, output_file, exclude_self)
