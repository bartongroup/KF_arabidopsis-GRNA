# KF_arabidopsis-GRNA
Scripts and commands for reproducing the analysis in "How well do RNA-Seq differential gene expression tools perform in a eukaryote with a complex transcriptome?", Froussios, Schurch, Mackinnon, Gierlinski, Duc, Simpson & Barton
doi: https://doi.org/10.1101/090753

# REQUIREMENTS:

How you install these requirements and make them accessible in your PATH is up to you.

1. Sun Grid Engine computing queue manager (some of our scripts internally submit tasks to a computing cluster)
2. R (3.2.2) with: lattice, ColourBrewer, PoissonSeq, edger, limma, deseq, deseq2, degseq, bayseq, ebseq, noiseq, samseq
3. PLplot
4. Perl (5.10.1) with: DRMAAc, PDL, PDL::Graphics::PLplot, DBI, Statistics::R, Devel::Size, Math::CDF, Time::HiRes
5. Python 2.7 with: numpy, scipy, drmaa, pysqlite
6. Python 3 with: pandas
7. STAR (2.5.0)
8. subread (1.6.1)
9. SAMtools (1.7)

Environment variables (assuming Unix/Bash):

	export PROJECTROOT=${HOME}/PROJECTS/AtGRNA
	export PERLROOT=${HOME}/perl5/perlbrew/perls/perl-5.10.1/
	export LOCALEXECROOT=${HOME}/local_installs

	export LD_LIBRARY_PATH=${SGE_ROOT}/lib/`${SGE_ROOT}/util/arch`
	export DRMAA_LIBRARY_PATH=${LD_LIBRARY_PATH}/libdrmaa.sa
	export PYTHONPATH=${PROJECTROOT}/analysis:${PROJECTROOT}/analysis/grnascripts/Modules
	export PERL5LIB=${PROJECTROOT}/analysis/grnascripts/Modules:${PROJECTROOT}/analysis/rnascripts:${PERLROOT}/lib/5.10.1:${PERL5ROOT}/lib/site_perl/5.10.1:${PERL5LIB}

You will likely need to change the paths in the above variables to match your system, especially the first three.

# PROCESSING STEPS

Commands are executed from the base directory of the code distribution. 

## Setup

First, create a project directory and navigate into it. Then clone this repository, get the Araport 11 *Arabidopsis thaliana* annotation and the TAIR10 *A. thaliana* genome assembly, and the raw read files from ArrayExpress.

	mkdir AtGRNA
	cd AtGRNA
	git clone git@github.com:bartongroup/KF_arabidopsis-GRNA.git ./
	mkdir ./genome
	ln -s ~/GENOMES/Arabidopsis_thaliana/Annotations/Araport-11/Araport11_genes.20151202.gff3 ./genome
	ln -s ~/GENOMES/Arabidopsis_thaliana/Annotations/Araport-11/Araport11_genes.20151202.gtf ./genome
	ln -s ~/GENOMES/Arabidopsis_thaliana/Annotations/Araport-11/At_tair10.fa ./genome
	mkdir fastq-unlev

Place or link the FASTQ files into the above subdirectory. Rename them AtWT_1.fastq, ...etc..., AtWT_17.fastq . Then we need to modify the chromosome labels in the annotation files to match those used in the FASTQ files.

	perl -e 'while(<>){~s/^Chr//g; print}' ./genome/Araport11_genes.20151202.gff3 > tmp.gff3
	perl -e 'while(<>){~s/^C/Pt/g; print}' tmp.gff3 > tmp2.gff3
	perl -e 'while(<>){~s/^M/Mt/g; print}' tmp2.gff3 > ./genome/At_araport11.gff3
	perl -e 'while(<>){~s/^Chr//g; print}' ./genome/Araport11_genes.20151202.gtf > tmp.gtf
	perl -e 'while(<>){~s/^C/Pt/g; print}' tmp.gtf > tmp2.gtf
	perl -e 'while(<>){~s/^M/Mt/g; print}' tmp2.gtf > ./genome/At_araport11.gtf
	rm ./tmp*

## STAR alignment and featureCounts gene expression

Make the STAR index, align the raw FASTQdata, count the reads mapping to genes.

	mkdir ./starindex
	mkdir ./starindex/At_araport11_99bp
	qsub -V -cwd -pe smp 8 -N staridx -b y STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./starindex/At_araport11_99bp --genomeFastaFiles ./genome/At_tair10.fa --sjdbGTFfile ./genome/At_araport11.gtf --sjdbOverhang 99
	mkdir sam-unlev
	qsub -V -cwd -pe smp 16 -N staral-unlev -t 1-17 ./analysis/star_array.sh -X STAR -c AtWT -i fastq-unlev -a At_araport11.gtf -n At_araport11_99bp -o sam-unlev -t 16
	python3 ./analysis/fileutilities.py T ./sam-unlev/ --dir final | python3 analysis/sequtilities.py P --StarFinalLogs tab -v > ./sam-unlev/AtWT-unlev_summary.tsv
	mkdir ./featurecounts-unlev
	> qsub -V cwd -pe smp 4 -N featCnt-unlev -t 1-17 ./analysis/featureCounts_array.sh -X featureCounts -c AtWT -a At_araport11.gtf -i sam-unlev -o featurecounts-unlev -f exon -g gene_id -r 2 -t 4

Collect all the read counts in a single table to use downstream

	mkdir ./combined_counts
	python3 ./analysis/fileutilities.py T ./featurecounts-unlev --dir "tsv$" | perl -e 'while(<>){~s/__Aligned.out.bam.exon\n/\n/;print}' > ./featurecounts-unlev/AtWT-unlev.list
	python3 ./analysis/fileutilities.py L ./featurecounts-unlev/AtWT-unlev.list --cols \-1 -li | perl -e 'while(<>){~s/_\|-1//g;print}' > combined_counts/AtWT-unlev_raw.tsv


## Calculating the correlations between replicates

	mkdir ./results
	Rscript ./analysis/correlations.R AtWT-unlev_raw.tsv > results/AtWT-unlev_raw.tsv_cors.csv
	Rscript ./analysis/levelplot.R results/AtWT-unlev_raw.tsv_cors.tsv

## FDR WT v WT Differential Gene expression (DGE)

This section uses code from the original [Schurch et. al. 2016](https://rnajournal.cshlp.org/content/22/6/839.long) paper to peform WT vs WT DGE calculations for selections of 3-7 replicates in each 'condition', boortstrapped 100 times.First we get the gene identifiers, setup some sub directories and the readcounts for each replicate.

	python3 ./analysis/fileutilities.py T ./combined_counts/AtWT-unlev_raw.tsv -r --cols 0 > At_genes-h.txt
	mkdir DGE_FDR_in
	mkdir ./DGE_FDR_in/AtWTa
	mkdir ./DGE_FDR_in/AtWTb
	mkdir DGE_FDR_out
	python3 ./analysis/fileutilities.py L ./featurecounts-unlev/AtWT-unlev.list -rl --cols 0 6 --out DGE_FDR_in/AtWTa "" '.gbgout'
	python3 ./analysis/fileutilities.py L ./featurecounts-unlev/AtWT-unlev.list -rl --cols 0 6 --out DGE_FDR_in/AtWTb "" '.gbgout'

Then we exclude the "bad" replicate 11.

	rm DGE_FDR_in/AtWTa/AtWT_11.gbgout DGE_FDR_in/AtWTb/AtWT_11.gbgout

Then we perform the DGE bootstrapping...

	bash ./analysis/loop.sh -f 2 -t 7 python ./analysis/grnascripts/Bootstrapping/generic_wrapper.py -o ${PROJECTROOT}/DGE_FDR_out/AtWT_edgeR_100_{val}.sqlite -l ${PROJECTROOT}/DGE_FDR_out/AtWT_edgeR_100_{val}.log -k {val} -r ${PROJECTROOT}/analysis/grnascripts/DE_tool_scripts/edgeR.R --tmpdir=${PROJECTROOT}/NOBACK/tmp-edgeR -n package:default -d ${PROJECTROOT}/DGE_FDR_in -a ${PROJECTROOT}/genome/At_araport11.gff3 -b 100 --precounts --gbgfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/group_by_gene.pl --agncfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/add_gene_name_column.pl --samtoolspath `which samtools` --Rpath `which Rscript`
	bash ./analysis/loop.sh -f 3 -t 7 python ./analysis/grnascripts/Bootstrapping/generic_wrapper.py -o ${PROJECTROOT}/DGE_FDR_out/AtWT_edgeRglm_100_{val}.sqlite -l ${PROJECTROOT}/DGE_FDR_out/AtWT_edgeRglm_100_{val}.log -k {val} -r ${PROJECTROOT}/analysis/grnascripts/DE_tool_scripts/edgeRglm.R --tmpdir=${PROJECTROOT}/NOBACK/tmp-edgeRglm -n package:default -d ${PROJECTROOT}/DGE_FDR_in -a ${PROJECTROOT}/genome/At_araport11.gff3 -b 100 --precounts --gbgfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/group_by_gene.pl --agncfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/add_gene_name_column.pl --samtoolspath `which samtools` --Rpath `which Rscript`
	bash ./analysis/loop.sh -f 3 -t 7 python ./analysis/grnascripts/Bootstrapping/generic_wrapper.py -o ${PROJECTROOT}/DGE_FDR_out/AtWT_deseq2_100_{val}.sqlite -l ${PROJECTROOT}/DGE_FDR_out/AtWT_deseq2_100_{val}.log -k {val} -r ${PROJECTROOT}/analysis/grnascripts/DE_tool_scripts/deseq2.R --tmpdir=${PROJECTROOT}/NOBACK/tmp-deseq2 -n package:default -d ${PROJECTROOT}/DGE_FDR_in -a ${PROJECTROOT}/genome/At_araport11.gff3 -b 100 --precounts --gbgfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/group_by_gene.pl --agncfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/add_gene_name_column.pl --samtoolspath `which samtools` --Rpath `which Rscript`
	bash ./analysis/loop.sh -f 3 -t 7 python ./analysis/grnascripts/Bootstrapping/generic_wrapper.py -o ${PROJECTROOT}/DGE_FDR_out/AtWT_deseq_100_{val}.sqlite -l ${PROJECTROOT}/DGE_FDR_out/AtWT_deseq_100_{val}.log -k {val} -r ${PROJECTROOT}/analysis/grnascripts/DE_tool_scripts/deseq.R --tmpdir=${PROJECTROOT}/NOBACK/tmp-deseq -n package:default -d ${PROJECTROOT}/DGE_FDR_in -a ${PROJECTROOT}/genome/At_araport11.gff3 -b 100 --precounts --gbgfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/group_by_gene.pl --agncfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/add_gene_name_column.pl --samtoolspath `which samtools` --Rpath `which Rscript`
	bash ./analysis/loop.sh -f 3 -t 7 python ./analysis/grnascripts/Bootstrapping/generic_wrapper.py -o ${PROJECTROOT}/DGE_FDR_out/AtWT_bayseq_100_{val}.sqlite -l ${PROJECTROOT}/DGE_FDR_out/AtWT_bayseq_100_{val}.log -k {val} -r ${PROJECTROOT}/analysis/grnascripts/DE_tool_scripts/bayseq.R --tmpdir=${PROJECTROOT}/NOBACK/tmp-bayseq -n package:default -d ${PROJECTROOT}/DGE_FDR_in -a ${PROJECTROOT}/genome/At_araport11.gff3 -b 100 --precounts --gbgfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/group_by_gene.pl --agncfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/add_gene_name_column.pl --samtoolspath `which samtools` --Rpath `which Rscript`
	bash ./analysis/loop.sh -f 3 -t 7 python ./analysis/grnascripts/Bootstrapping/generic_wrapper.py -o ${PROJECTROOT}/DGE_FDR_out/AtWT_degseq_100_{val}.sqlite -l ${PROJECTROOT}/DGE_FDR_out/AtWT_degseq_100_{val}.log -k {val} -r ${PROJECTROOT}/analysis/grnascripts/DE_tool_scripts/degseq.R --tmpdir=${PROJECTROOT}/NOBACK/tmp-degseq -n package:default -d ${PROJECTROOT}/DGE_FDR_in -a ${PROJECTROOT}/genome/At_araport11.gff3 -b 100  --precounts --gbgfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/group_by_gene.pl --agncfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/add_gene_name_column.pl --samtoolspath `which samtools` --Rpath `which Rscript`
	bash ./analysis/loop.sh -f 3 -t 7 python ./analysis/grnascripts/Bootstrapping/generic_wrapper.py -o ${PROJECTROOT}/DGE_FDR_out/AtWT_ebseq_100_{val}.sqlite -l ${PROJECTROOT}/DGE_FDR_out/AtWT_ebseq_100_{val}.log -k {val} -r ${PROJECTROOT}/analysis/grnascripts/DE_tool_scripts/ebseq.R --tmpdir=${PROJECTROOT}/NOBACK/tmp-ebseq -n package:default -d ${PROJECTROOT}/DGE_FDR_in -a ${PROJECTROOT}/genome/At_araport11.gff3 -b 100 --precounts --gbgfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/group_by_gene.pl --agncfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/add_gene_name_column.pl --samtoolspath `which samtools` --Rpath `which Rscript`
	bash ./analysis/loop.sh -f 3 -t 7 python ./analysis/grnascripts/Bootstrapping/generic_wrapper.py -o ${PROJECTROOT}/DGE_FDR_out/AtWT_limma_100_{val}.sqlite -l ${PROJECTROOT}/DGE_FDR_out/AtWT_limma_100_{val}.log -k {val} -r ${PROJECTROOT}/analysis/grnascripts/DE_tool_scripts/limma.R --tmpdir=${PROJECTROOT}/NOBACK/tmp-limma -n package:default -d ${PROJECTROOT}/DGE_FDR_in -a ${PROJECTROOT}/genome/At_araport11.gff3 -b 100 --precounts --gbgfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/group_by_gene.pl --agncfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/add_gene_name_column.pl --samtoolspath `which samtools` --Rpath `which Rscript`
	bash ./analysis/loop.sh -f 3 -t 7 python ./analysis/grnascripts/Bootstrapping/generic_wrapper.py -o ${PROJECTROOT}/DGE_FDR_out/AtWT_samseq_100_{val}.sqlite -l ${PROJECTROOT}/DGE_FDR_out/AtWT_samseq_100_{val}.log -k {val} -r ${PROJECTROOT}/analysis/grnascripts/DE_tool_scripts/samseq.R --tmpdir=${PROJECTROOT}/NOBACK/tmp-samseq -n package:default -d ${PROJECTROOT}/DGE_FDR_in -a ${PROJECTROOT}/genome/At_araport11.gff3 -b 100 --precounts --gbgfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/group_by_gene.pl --agncfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/add_gene_name_column.pl --samtoolspath `which samtools` --Rpath `which Rscript`
	bash ./analysis/loop.sh -f 3 -t 7 python ./analysis/grnascripts/Bootstrapping/generic_wrapper.py -o ${PROJECTROOT}/DGE_FDR_out/AtWT_poissonseq_100_{val}.sqlite -l ${PROJECTROOT}/DGE_FDR_out/AtWT_poissonseq_100_{val}.log -k {val} -r ${PROJECTROOT}/analysis/grnascripts/DE_tool_scripts/poissonseq.R --tmpdir=${PROJECTROOT}/NOBACK/tmp-poissonseq -n package:default -d ${PROJECTROOT}/DGE_FDR_in -a ${PROJECTROOT}/genome/At_araport11.gff3 -b 100 --precounts --gbgfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/group_by_gene.pl --agncfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/add_gene_name_column.pl --samtoolspath `which samtools` --Rpath `which Rscript`
	bash ./analysis/loop.sh -f 3 -t 7 python ./analysis/grnascripts/Bootstrapping/generic_wrapper.py -o ${PROJECTROOT}/DGE_FDR_out/AtWT_noiseq_100_{val}.sqlite -l ${PROJECTROOT}/DGE_FDR_out/AtWT_noiseq_100_{val}.log -k {val} -r ${PROJECTROOT}/analysis/grnascripts/DE_tool_scripts/noiseq.R --tmpdir=${PROJECTROOT}/NOBACK/tmp-noiseq -n package:default -d ${PROJECTROOT}/DGE_FDR_in -a ${PROJECTROOT}/genome/At_araport11.gff3 -b 100 --precounts --gbgfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/group_by_gene.pl --agncfile ${PROJECTROOT}/analysis/grnascripts/Bootstrapping/add_gene_name_column.pl --samtoolspath `which samtools` --Rpath `which Rscript`

With that done we can compute the FDR stats and plot them for the paper...

	mkdir ./DGE_FDR_powerstats
	source ./analyis/loop.sh -v bayseq -v deseq1 -v deseq2 -v degseq -v edger -v edgerglm -v ebseq -v limma -v noiseq -v poissonseq -v samseq "python ./analysis/grnascripts/DE_tool_comparison/make_powerstats_db.pl -test={val} -testinfofile=${PROJECTROOT}/analysis/DGE_FDR_tests_samecond.txt -script=${PROJECTROOT}/analysis/grnascripts/DE_tool_comparison/one_bs_powerstats.pl -genlist=At_genes-h.txt -powerdir=DGE_FDR_powerstats -reffcfile=${PROJECTROOT}/analysis/DGE_FDR_dummy-truth.tsv -maxn=7"
	perl ./analysis/grnascripts/Plotting/compare_null.pl -psfile=results/dge-fdr_bay-deg-de1-de2-eb-edger-eglm-limma-noi-pois-sam.ps -testinfofile=${PROJECTROOT}/analysis/DGE_FDR_tests_samecond.txt -dir=DGE_FDR_powerstats -maxy=1 -tests=bayseq,degseq,deseq1,deseq2,ebseq,edger,edgerglm,limma,noiseq,poissonseq,samseq

## Normalising replicates by sampling reads

Here we downsample the replicates to the lowest replicate so that we can use statistical tests of goodness-of-fit to distributions that that require integer counts. We then Align the data, and get the gene counts (cluding subsetting the replicates for the 'less noisy' set comparison we make in the paper).

	mkdir fastq-lev
	> python3 ./analysis/fileutilities.py T ./sam-unlev/AtWT-unlev_summary.tsv --cols 0 2 -rl | perl -e 'while(<>){~s/_/\t/g;print}' > ./combined_counts/AtWT_readcount.dat
	> qsub -cwd -V -t 1-17 ./analysis/levelling_array.sh -X analysis/levelling.pl -c AtWT -p 1 -i fastq-unlev -o fastq-lev -n combined_counts/AtWT_readcount.dat
	mkdir sam-lev
	qsub -V -cwd -pe smp 16 -N staral-lev -t 1-17 ./analysis/star_array.sh -X STAR -c AtWT -i fastq-lev -a At_araport11.gtf -n At_araport11_99bp -o sam-lev -t 16
	mkdir featurecounts-lev
	qsub -V cwd -pe smp 4 -N featCnt-lev -t 1-17 ./analysis/featureCounts_array.sh -X featureCounts -c AtWT -a At_araport11.gtf -i sam-lev -o featurecounts-lev -f exon -g gene_id -r 2 -t 4
	python3 ./analysis/fileutilities.py T ./featurecounts-lev --dir "tsv$" | python3 analysis/fileutilities.py P --cols \-1 -lri > ./combined_counts/AtWT-lev_raw.tsv

Collect the counts of the less "noisy" replicates, and for an equally-sized "control" group that includes the "noisy" replicates, to exclude the influence of the number of replicates.

	python3 ./analysis/fileutilities.py T ./combined_counts/AtWT-lev_raw.tsv -ri --cols 1 2 3 4 5 6 7 15 16 17 > ./combined_counts/AtWT-levN_raw.tsv
	python3 ./analysis/fileutilities.py T ./combined_counts/AtWT-lev_raw.tsv -ri --cols 8 9 10 11 12 13 14 15 16 17 > ./combined_counts/AtWT-levC_raw.tsv

## Goodness-of-fit distribution tests

On this data we not run the good-ness of fit tests for the distributions and plot the results for the figures in the paper.

	ln -s /.analysis/defs.dat ./
	qsub -V -cwd -b y perl ./analysis/rnascripts/one_nb_test.pl -batch=1 -batchsize=33851 -outfile=AtWT-lev_nbtest_m_deseq.dat -stat=m -cond=AtWT-lev -norm=deseq -nonzero=1 -type=raw -test=meintanis -ncrit=30 -maxn=10000000
	bash ./analysis/loop.sh -v nb -v norm -v lnorm -v pois 'perl ./analysis/rnascripts/distribution_test.pl -dist={val} -psfile=results/AtWT_{val}_test.ps -cond=AtWT-lev -multicor=bh'

Do the same for the "non-noisy" and "control" subsets as well

	qsub -V -cwd -b y perl ./analysis/rnascripts/one_nb_test.pl -batch=1 -batchsize=33851 -outfile=AtWT-levN_nbtest_m_deseq.dat -stat=m -cond=AtWT-levN -norm=deseq -nonzero=1 -type=raw -test=meintanis -ncrit=30 -maxn=10000000
	bash ./analysis/loop.sh -v nb -v norm -v lnorm -v pois 'perl ./analysis/rnascripts/distribution_test.pl -dist={val} -psfile=results/AtWTn_{val}_test.ps -cond=AtWT-levN -multicor=bh'
	qsub -V -cwd -b y perl ./analysis/rnascripts/one_nb_test.pl -batch=1 -batchsize=33851 -outfile=AtWT-levC_nbtest_m_deseq.dat -stat=m -cond=AtWT-levC -norm=deseq -nonzero=1 -type=raw -test=meintanis -ncrit=30 -maxn=10000000
	bash ./analysis/loop.sh -v nb -v norm -v lnorm -v pois 'perl ./analysis/rnascripts/distribution_test.pl -dist={val} -psfile=results/AtWTc_{val}_test.ps -cond=AtWT-levC -multicor=bh'

And thats it.
