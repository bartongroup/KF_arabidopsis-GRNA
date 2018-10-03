source ~/local_installs/miniconda3/bin/activate atgrna
qsub -cwd -V -R y -pe smp 4 -N featCnt-lev -t 1-17 ./analysis/featureCounts_array.sh -X featureCounts -c AtWT -a At_araport11.gtf -i sam-lev -o featurecounts-lev -f exon -g gene_id -r 2 -t 4
