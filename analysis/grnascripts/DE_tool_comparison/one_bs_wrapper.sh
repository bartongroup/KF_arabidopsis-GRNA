#!/bin/bash

export PROJECTROOT=${HOME}/PROJECTS/AtGRNA
export PERLROOT=${HOME}/perl5/perlbrew/perls/perl-5.10.1/
export LOCALEXECROOT=${HOME}/local_installs

cd ${PROJECTROOT}
script=${PROJECTROOT}/analysis/grnascripts/DE_tool_comparison/one_bs_powerstats.pl

source ${HOME}/perl5/perlbrew/etc/bashrc
perlbrew use 5.10.1

export LD_LIBRARY_PATH=${SGE_ROOT}/lib/`${SGE_ROOT}/util/arch`
export DRMAA_LIBRARY_PATH=${LD_LIBRARY_PATH}/libdrmaa.so
export PYTHONPATH=${PROJECTROOT}/analysis:${PROJECTROOT}/analysis/grnascripts/Modules
export PERL5LIB=${PROJECTROOT}/analysis/grnascripts/Modules;${PERLROOT}/lib/5.10.1;${PERLROOT}/lib/site_perl/5.10.1

source ~/local_installs/miniconda3/bin/activate atgrna

command="perl ${script} $@"
${command}

#perl ./analysis/grnascripts/DE_tool_comparison/one_bs_powerstats.pl
