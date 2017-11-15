#!/bin/bash

# Send array job for levelling.pl. 
# SGE_TASK_ID <-> replicate number

function usage() {
    echo "Usage:"
    echo "      $0 -c PREFIX -p PAIRED [-i FASTQDIR] [-o OUTDIR] [-n COUNTSFILE] [-X path/levelling.pl] [-D PROJECTDIR] [-L]"
    echo "All directories must be relative to PROJECTDIR."
    exit 1
}

# Set defaults.
project=$(pwd)
tool='./levelling.pl'
cond=''
rep=''
countsfile='combined_counts/readcount.dat'
fastqdir='fastq-unlev'
leveldir='fastq-lev'
pair=0
log=false
scriptsdir='./'


# Parse options.
while getopts 'c:r:n:i:o:p:X:D:L' flag; do
  case "${flag}" in
    c) cond="${OPTARG}" ;;
    i) fastqdir="${OPTARG}" ;;
    o) leveldir="${OPTARG}" ;;
    n) countsfile="${OPTARG}" ;;
    p) pair="${OPTARG}" ;;
    X) tool="${OPTARG}" ;;
    D) project="${OPTARG}" ;;
    L) log=true ;;
    *) usage ;;
  esac
done

# No point going further without a prefix.
[[ -z "$cond" ]] && usage


# Log and Execute

command="$tool -c $cond -r $SGE_TASK_ID -p $pair -i ${project}/${fastqdir} -o ${project}/${leveldir} -n ${project}/${countsfile}"
if [ "$log" = true ] ; then  
    python "${scriptsdir}/mylogs.py" "$command"
fi

$command
