if [ "$#" -lt 2 ]
then
    echo "Usage: $0 <url and destination directory>"
    exit 1
fi

url=$1
destination_directory=$2
uncompress=$3
echo "Downloading the sequencing data files..."
wget -P ~/decont/${destination_directory} $url

if [ "$uncompress" = "yes" ]
then
    wget -P ~/decont/${destination_directory} https://bioinformatics.cnio.es/data/courses/decont/contaminants.fasta.gz | gunzip ~/decont/${destination_directory}/contaminants.fasta.gz
    awk '/small nuclear RNA/ {flag=1; next} /^>/ {flag=0} !flag' ~/decont/${destination_directory}/contaminants.fasta > ~/decont/${destination_directory}/filtered_contaminants.fasta
fi


# This script should download the file specified in the first argument ($1),
# place it in the directory specified in the second argument ($2),
# and *optionally*:
# - uncompress the downloaded file with gunzip if the third
#   argument ($3) contains the word "yes"
# - filter the sequences based on a word contained in their header lines:
#   sequences containing the specified word in their header should be **excluded**
#
# Example of the desired filtering:
#
#   > this is my sequence
#   CACTATGGGAGGACATTATAC
#   > this is my second sequence
#   CACTATGGGAGGGAGAGGAGA
#   > this is another sequence
#   CCAGGATTTACAGACTTTAAA
#
#   If $4 == "another" only the **first two sequence** should be output
