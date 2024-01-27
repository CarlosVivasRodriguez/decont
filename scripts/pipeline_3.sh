#!/bin/bash

# Download all the files specified in data/filenames
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
	# Obtén el nombre del archivo sin la ruta ni la extensión
	filename=$(basename "$filepath")
	# Utiliza sed para realizar la transformación del nombre del archivo
	sample_id=$(echo "$filename" | sed -E 's/-12\.5dpp\.1\.[12]s_sRNA.*//;s/\..*//')
	# Ejecuta el script de fusión
	bash scripts/merge_fastqs.sh data out/merged "$sample_id"
done

# Ejecutar cutadapt para todos los archivos fusionados
mkdir -p log/cutadapt  # Crear la carpeta cutadapt dentro de log
for merged_file in out/merged/*.fastq.gz; do
    trimmed_file="out/trimmed/$(basename "$merged_file" .fastq.gz).trimmed.fastq.gz"
    cutadapt -m 18 -a TGGAATTCTCGGGTGCCAAGG --discard-untrimmed \
        -o "$trimmed_file" "$merged_file" > "log/cutadapt/$(basename "$merged_file" .fastq.gz).log"

    # Añadir información relevante al archivo de registro principal
    echo "Cutadapt results for $(basename "$merged_file"):" >> "$log_file"
    grep -E 'Reads with adapters|Total basepairs' "log/cutadapt/$(basename "$merged_file" .fastq.gz).log" >> "$log_file"
done


# Ejecutar STAR para todos los archivos recortados
mkdir -p out/star
for trimmed_file in out/trimmed/*.fastq.gz; do
    sid=$(basename "$trimmed_file" .trimmed.fastq.gz)
    output_directory="out/star/$sid"
    mkdir -p "$output_directory"
    STAR --runThreadN 4 --genomeDir res/contaminants_idx \
         --outReadsUnmapped Fastx --readFilesIn "$trimmed_file" \
         --readFilesCommand gunzip -c --outFileNamePrefix "$output_directory/" > temp_star.log

    # Añadir información relevante al archivo de registro
    echo "STAR results for $sid:" >> "$log_file"
    grep -E 'Uniquely mapped reads %|Number of reads mapped to multiple loci|% of reads mapped to too many loci' temp_star.log >> "$log_file"
done

# Limpiar archivos temporales de log
rm temp_cutadapt.log temp_star.log
