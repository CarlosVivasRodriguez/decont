#!/bin/bash

# Download all the files specified in data/urls
for url in $(<data/urls) #TODO
do
    bash scripts/download.sh "$url" data
done

# Download the contaminants fasta file, uncompress it, and
# filter to remove all small nuclear RNAs
contaminants_url="https://bioinformatics.cnio.es/data/courses/decont/contaminants.fasta.gz"
bash scripts/download.sh "$contaminants_url" res yes #TODO

# Index the contaminants file
bash scripts/index.sh res/contaminants.fasta res/contaminants_idx

# Merge the samples into a single file
for filepath in data/*fastq.gz; do
    # Obtain the filename without path and extension
    filename=$(basename "$filepath")
    # Use sed to transform the filename
    sample_id=$(echo "$filename" | sed -E 's/-12\.5dpp\.1\.[12]s_sRNA.*//;s/\..*//')
    # Execute the merge script
    bash scripts/merge_fastqs.sh data out/merged "$sample_id"
done

# Execute cutadapt for all merged files
log_file="log/pipeline.log"
mkdir -p log/cutadapt  # Create log/cutadapt directory if it doesn't exist
mkdir -p out/trimmed
for merged_file in out/merged/*.fastq.gz; do
    trimmed_file="out/trimmed/$(basename "$merged_file" .fastq.gz).trimmed.fastq.gz"
    cutadapt -m 18 -a TGGAATTCTCGGGTGCCAAGG --discard-untrimmed \
        -o "$trimmed_file" "$merged_file" > log/cutadapt/$(basename "$merged_file" .fastq.gz).log

    # Add relevant information to the log file
    echo "Cutadapt results for $(basename "$merged_file"):" >> "$log_file"
    grep -E 'Reads with adapters|Total basepairs' log/cutadapt/$(basename "$merged_file" .fastq.gz).log >> "$log_file"
done

# Execute STAR for all trimmed files
mkdir -p out/star
for trimmed_file in out/trimmed/*.fastq.gz; do
    sid=$(basename "$trimmed_file" .trimmed.fastq.gz)
    output_directory="out/star/$sid"
    mkdir -p "$output_directory"
    STAR --runThreadN 4 --genomeDir res/contaminants_idx \
         --outReadsUnmapped Fastx --readFilesIn "$trimmed_file" \
         --readFilesCommand gunzip -c --outFileNamePrefix "$output_directory/" > temp_star.log

    # Add relevant information to the log file
    echo "STAR results for $sid:" >> "$log_file"
    grep -E 'Uniquely mapped reads %|Number of reads mapped to multiple loci|% of reads mapped to too many loci' temp_star.log >> "$log_file"
done

# Clean up temporary log files
rm temp_star.log

