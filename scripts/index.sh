#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <genome_file> <output_directory>"
    exit 1
fi

genome_file=$1
output_directory=$2

mkdir -p "$output_directory"

STAR --runThreadN 4 --runMode genomeGenerate \
     --genomeDir "$output_directory" \
     --genomeFastaFiles "$genome_file" \
     --genomeSAindexNbases 9
