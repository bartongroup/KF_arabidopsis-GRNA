source ~/local_installs/miniconda3/bin/activate atgrna
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./starindex/At_araport11_99bp --genomeFastaFiles ./genome/At_tair10.fa --sjdbGTFfile ./genome/At_araport11.gtf --sjdbOverhang 99
