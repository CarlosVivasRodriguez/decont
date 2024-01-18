#!/bin/bash

# This script should index the genome file specified in the first argument ($1),
# creating the index in a directory specified by the second argument ($2).

# The STAR command is provided for you. You should replace the parts surrounded
# by "<>" and uncomment it.

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <genome_file> <output_directory>"
    exit 1
fi

# Assign the input arguments to variables
genome_file=$1 
output_directory=$2 

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Run STAR command for genome indexing
STAR --runThreadN 4 --runMode genomeGenerate \
     --genomeDir "$output_directory" \
     --genomeFastaFiles "$genome_file" \
     --genomeSAindexNbases 9
