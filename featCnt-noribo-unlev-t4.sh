source ~/local_installs/miniconda3/bin/activate atgrna
source ./analysis/featureCounts_array.sh -X featureCounts -c AtWT -a At_araport11.gtf -i sam-noribo-unlev -o featurecounts-noribo-unlev -f exon -g gene_id -r 2 -t 4
