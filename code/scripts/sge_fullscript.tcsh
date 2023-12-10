#!/bin/tcsh

##
## This script is for submitting fullscript.sh jobs to 
## the sge cluster in the physics node computing system
## 
## Usage:
##   user> qsub sge_fullscript.sh inputfile.in
##

#$ -S /bin/tcsh
#$ -j y
#$ -q all.q

source ~/.cshrc

echo "################################################"
env
echo "################################################"

${MG_HOME}/MIW_codebase/scripts/fullscript.sh ${1}
