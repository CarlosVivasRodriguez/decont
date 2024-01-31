#!/bin/bash

for url in $(<data/urls)
do
    echo "Downloading: $url"
    bash scripts/download.sh "$url" data
done

contaminants_url="https://bioinformatics.cnio.es/data/courses/decont/contaminants.fasta.gz"
echo "Downloading contaminants file: $contaminants_url"
bash scripts/download.sh "$contaminants_url" res yes

bash scripts/index.sh res/contaminants.fasta res/contaminants_idx

list_of_sample_ids=$(ls -1 data/*.fastq.gz | grep -E "/[A-Z][^\]*$" | xargs -I {} basename {} | cut -d"-" -f1| sort | uniq)

for sid in $list_of_sample_ids; do
    bash scripts/merge_fastqs.sh data out/merged "$sid"
done

log_file="log/pipeline.log"
mkdir -p log/cutadapt
mkdir -p out/trimmed
for merged_file in out/merged/*.fastq.gz; do
    trimmed_file="out/trimmed/$(basename "$merged_file" .fastq.gz).trimmed.fastq.gz"
    cutadapt -m 18 -a TGGAATTCTCGGGTGCCAAGG --discard-untrimmed \
        -o "$trimmed_file" "$merged_file" > log/cutadapt/$(basename "$merged_file" .fastq.gz).log


    echo "Cutadapt results for $(basename "$merged_file"):" >> "$log_file"
    grep -E 'Reads with adapters|Total basepairs' log/cutadapt/$(basename "$merged_file" .fastq.gz).log >> "$log_file"
done

mkdir -p out/star
for trimmed_file in out/trimmed/*.fastq.gz; do
    sid=$(basename "$trimmed_file" .trimmed.fastq.gz)
    output_directory="out/star/$sid"
    mkdir -p "$output_directory"
    STAR --runThreadN 4 --genomeDir res/contaminants_idx \
         --outReadsUnmapped Fastx --readFilesIn "$trimmed_file" \
         --readFilesCommand gunzip -c --outFileNamePrefix "$output_directory/" > temp_star.log

    echo "STAR results for $sid:" >> "$log_file"
    grep -E 'Uniquely mapped reads %|Number of reads mapped to multiple loci|% of reads mapped to too many loci' temp_star.log >> "$log_file"
done

rm temp_star.log
