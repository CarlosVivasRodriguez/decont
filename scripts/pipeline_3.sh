# Descargar todos los archivos especificados en data/urls
for url in $(<data/urls)
do
    bash scripts/download.sh "$url" data
done


# Descargar el archivo fasta de contaminantes, descomprimirlo y
# filtrar para eliminar todos los ARN nucleares pequeños
contaminants_url="https://bioinformatics.cnio.es/data/courses/decont/contaminants.fasta.gz"
bash scripts/download.sh "$contaminants_url" res yes
# TODO: Agregar código para filtrar ARN nucleares pequeños si es necesario

# Indexar el archivo de contaminantes
bash scripts/index.sh res/contaminants.fasta res/contaminants_idx

# Fusionar todas las muestras en un solo archivo
for fastq_file in data/*.fastq.gz; do
    sid=$(basename "$fastq_file" .fastq.gz)
    bash scripts/merge_fastqs.sh data out/merged "$sid"
done

# Ejecutar cutadapt para todos los archivos fusionados
log_file="log/pipeline.log"
mkdir -p out/trimmed
for merged_file in out/merged/*.fastq.gz; do
    trimmed_file="out/trimmed/$(basename "$merged_file" .fastq.gz).trimmed.fastq.gz"
    cutadapt -m 18 -a TGGAATTCTCGGGTGCCAAGG --discard-untrimmed \
        -o "$trimmed_file" "$merged_file" > temp_cutadapt.log

    # Añadir información relevante al archivo de registro
    echo "Cutadapt results for $(basename "$merged_file"):" >> "$log_file"
    grep -E 'Reads with adapters|Total basepairs' temp_cutadapt.log >> "$log_file"
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
