source ~/local_installs/miniconda3/bin/activate atgrna
source ./analysis/star_array.sh -X STAR -c AtWT -i fastq-lev -a At_araport11.gtf -n At_araport11_99bp -o sam-lev -t 16
