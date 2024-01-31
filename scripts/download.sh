if [ "$#" -lt 2 ]
then
    echo "Usage: $0 <url and destination directory>"
    exit 1
fi
url=$1
destination_directory=$2
uncompress=$3
wget -P ~/decont/${destination_directory} $url
if [ "$uncompress" = "yes" ]
then
    wget -P ~/decont/${destination_directory} https://bioinformatics.cnio.es/data/courses/decont/contaminants.fasta.gz
    gunzip ~/decont/${destination_directory}/contaminants.fasta.gz | awk '/small nuclear RNA/ {flag=1; next} /^>/ {flag=0} !flag' ~/decont/${destination_directory}/contaminants.fasta > ~/decont/${destination_directory}/contaminants.fasta
fi
