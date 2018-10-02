source ~/perl5/perlbrew/etc/bashrc
perlbrew use 5.10.1

export PROJECTROOT=${HOME}/PROJECTS/AtGRNA
export PERLROOT=${HOME}/perl5/perlbrew/perls/perl-5.10.1/
export LOCALEXECROOT=${HOME}/local_installs

export LD_LIBRARY_PATH=${SGE_ROOT}/lib/`${SGE_ROOT}/util/arch`
export DRMAA_LIBRARY_PATH=${LD_LIBRARY_PATH}/libdrmaa.so
export PYTHONPATH=${PROJECTROOT}/analysis:${PROJECTROOT}/analysis/grnascripts/Modules
export PERL5LIB=${PROJECTROOT}/analysis/grnascripts/Modules:${PROJECTROOT}/analysis/rnascripts:${PERLROOT}/lib/5.10.1:${PERLROOT}/lib/site_perl/5.10.1:${PERL5LIB}
