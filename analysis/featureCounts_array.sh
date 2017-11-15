#!/bin/bash

# Send array job for STAR. 
# SGE_TASK_ID <-> replicate number


# Set defaults.
log=false
threads=1
tool=/sw/opt/subread-1.4.6-p4-Linux-x86_64/bin/featureCounts
project=$(pwd)
feature='exon'
groupby='gene_id'
stranded=2
countsdir='raw_counts'
prefix=''
suffix='_Aligned.out.bam'
annotdir='genome'
annot=''
samdir='sam'
scriptsdir="./"

function usage() {
    echo "Usage:"
    echo "      $0 -c PREFIX -a ANNOTATION [-s SUFFIX] [-A ANNOTDIR] [-i SAMDIR] [-o COUNTSDIR] [-f FEATURE] [-g GROUPBY] [-r STRANDED] [-t THREADS] [-X path/featureCounts] [-D PROJECTDIR] [-L]"
    echo "Defaults:"
    echo "      -s $suffix -A $annotdir -i $samdir -o $countsdir -f $feature -g $groupby -r $stranded -t $threads -X $tool -D $project -L $log"
    echo "All directories must be relative to PROJECTDIR. -L requires ${scriptsdir}/mylogs.py"
    exit 1
}


# Parse options.
while getopts 'c:a:s:A:i:o:f:g:r:t:X:D:L' flag; do
  case "${flag}" in
	c) prefix="${OPTARG}" ;;
    a) annot="${OPTARG}" ;;
	s) suffix="${OPTARG}" ;;
    A) annotdir="${OPTARG}" ;;
    i) samdir="${OPTARG}" ;;
    o) countsdir="${OPTARG}" ;;
    f) feature="${OPTARG}" ;;
    g) groupby="${OPTARG}" ;;
    r) stranded="${OPTARG}" ;;
    t) threads="${OPTARG}" ;;
    X) tool="${OPTARG}" ;;
    D) project="${OPTARG}" ;;
    L) log=true ;;
    *) usage ;;
  esac
done

# No point going further without at least these.
[[ -z "$prefix" ]] && usage
[[ -z "$annot" ]] && usage

# Tidy up the arguments.
myparams="-a ${project}/${annotdir}/${annot} -o ${project}/${countsdir}/${prefix}_${SGE_TASK_ID}_${suffix}.${feature}.tsv -t $feature -g $groupby -s $stranded -T $threads"
miscparams="-p -P"
input="${project}/${samdir}/${prefix}_${SGE_TASK_ID}${suffix}"

# Log and Execute.

command="$tool $myparams $miscparams $input"
if [ "$log" = true ] ; then  
	python "${scriptsdir}/mylogs.py" "$command"
fi

$command
 
 
 
 
 

  
 

