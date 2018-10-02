#!/bin/bash

# Send array job for levelling.pl. 
# SGE_TASK_ID <-> replicate number

function usage() {
    echo "Usage:"
    echo "      $0 -c PREFIX -p PAIRED [-i FASTQDIR] [-o OUTDIR] [-n COUNTSFILE] [-X path/levelling.pl] [-D PROJECTDIR]"
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


# Parse options.
while getopts 'c:r:n:i:o:p:X:D' flag; do
  case "${flag}" in
    c) cond="${OPTARG}" ;;
    i) fastqdir="${OPTARG}" ;;
    o) leveldir="${OPTARG}" ;;
    n) countsfile="${OPTARG}" ;;
    p) pair="${OPTARG}" ;;
    X) tool="${OPTARG}" ;;
    D) project="${OPTARG}" ;;
    *) usage ;;
  esac
done

# No point going further without a prefix.
[[ -z "$cond" ]] && usage


# Execute
command="perl $tool -c $cond -r $SGE_TASK_ID -p $pair -i ${project}/${fastqdir} -o ${project}/${leveldir} -n ${project}/${countsfile}"

$command


#eof
