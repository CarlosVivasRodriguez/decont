if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <samples_directory> <output_directory> <sample_id>"
    exit 1
fi

samples_directory=$1
output_directory=$2
sample_id=$3

mkdir -p $output_directory

cat $samples_directory/${sample_id}*.fastq.gz > $output_directory/${sample_id}.fastq.gz
