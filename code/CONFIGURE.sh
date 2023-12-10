#! /bin/bash

echo ""
echo " Checking environment is compatible and ready for use..."
echo " "

## Python ####################################################################################

rm -f printversion.py scratch.txt
echo "#! /usr/bin/env python" > printversion.py
echo " " >> printversion.py
echo "import sys" >> printversion.py
echo "print(\"MIW-TAG %s\" % sys.version)" >> printversion.py
python printversion.py | sed "s/\(MIW-TAG [0-9].[0-9].[0-9]\).*/\1/g" > scratch.txt
PYTHON_VERSION=`grep -e "MIW-TAG" scratch.txt | sed "s/MIW-TAG //g" `
PYTHON_VERSION_PRIMARY=`echo ${PYTHON_VERSION} | sed "s/\([0-9]\).[0-9].[0-9]/\1/g"`
PYTHON_VERSION_SECONDARY=`echo ${PYTHON_VERSION} | sed "s/[0-9].\([0-9]\).[0-9]/\1/g"`
PYTHON_VERSION_TERTIARY=`echo ${PYTHON_VERSION} | sed "s/[0-9].[0-9].\([0-9]\)/\1/g"`
rm -r printversion.py scratch.txt

echo "PYTHON_VERSION: ${PYTHON_VERSION}"
if [ ${PYTHON_VERSION_PRIMARY} != 2 ] || [ ${PYTHON_VERSION_SECONDARY} -lt 6 ] ; then
  echo " ######## ERROR ######## "
  echo "MadGraph5 works only with python 2.6 or later (but not python 3.X)."
  echo "               Please upgrade your version of python.               "
  exit
fi

## MadGraph5 ####################################################################################

if [ -z ${MG_HOME} ] ; then
  echo "In your .cshrc file, please set MG_HOME as the location of the MadGraph5 directory."
  exit
else
  echo "MG_HOME:        ${MG_HOME}"
  if [ ! -d ${MG_HOME} ] ; then 
    echo "${MG_HOME} does not exist. Please create the directory before continuin."
    exit
  fi
fi

cd ${MG_HOME}
if [ ! -d MIW_codebase ] ; then 
  echo " ######## ERROR ######## "
  echo "The MIW_codebase directory needs to be placed within the MG_HOME directory."
  exit
else
  echo "                    --> MIW_codebase is in the correct location."
fi

if [ ! -e ${MG_HOME}/bin/mg5 ] ; then
  echo " ######## ERROR ######## "
  echo "The MadGraph5 executable 'mg5' should exist in ${MG_HOME}/bin/"
  exit
else
  echo "                    --> ${MG_HOME}/bin/mg5 exists."
fi

if [ ! -d ${MG_HOME}/pythia-pgs ] ; then
  echo " ######## ERROR ######## "
  echo "The pythia-pgs directory needs to be installed within the MG_HOME directory."
  exit
else
  echo "                    --> ${MG_HOME}/pythia-pgs exists."
fi


## Data Storage ####################################################################################

if [ -z ${DATA_HOME} ] ; then
  echo " ######## ERROR ######## "
  echo "In your .cshrc file, please set DATA_HOME as the directory you wish MadGraph5 data to be placed."
  echo "A default value of DATA_HOME=${MG_HOME} can be used."
  exit
else
  echo "DATA_HOME:      ${DATA_HOME}"
  if [ ! -d ${DATA_HOME} ] ; then 
    echo "${DATA_HOME} does not exist. Please create the directory before continuing."
    exit
  fi
fi
mkdir -p ${DATA_HOME}/Storage
mkdir -p ${DATA_HOME}/Processes
mkdir -p ${DATA_HOME}/Processes/sew_lhe_files

## CERN's ROOT package ####################################################################################

if [ -z ${ROOTSYS} ] ; then  
  echo " ######## ERROR ######## "
  echo "In your .cshrc file, please set ROOTSYS as the location of CERN's ROOT package."
  exit
else
  echo "ROOTSYS:        ${ROOTSYS}"
  if [ ! -f ${ROOTSYS} ] ; then 
    echo "${ROOTSYS} does not exist."
    exit
  fi
fi

rm -f scratch.txt
readelf -h ${ROOTSYS}/bin/root > scratch.txt
ELF_VERSION=`grep -e "Class:" scratch.txt | sed "s/.*\(ELF[0-9]*\)/\1/g"`
rm -f scratch.txt
if [ "${ELF_VERSION}" != "ELF64" ] ; then
  echo " ######## WARNING ######## "
  echo "It is highly recommended that a 64-bit architecture is used for CERN's ROOT package."
fi
echo "                    --> ${ELF_VERSION}"

## SUSY-HIT ####################################################################################

if [ -z ${SH_HOME} ] ; then
  echo " ######## ERROR ######## "
  echo "In your .cshrc file, please set SH_HOME as the location of the SUSY-HIT directory."
  exit
else
  echo "SH_HOME:        ${SH_HOME}"
  if [ ! -d ${SH_HOME} ] ; then 
    echo "${SH_HOME} does not exist. Please create the directory before continuing."
    exit
  fi
fi
if [ ! -e ${SH_HOME}/run ] ; then
  echo " ######## ERROR ######## "
  echo "The SUSY-HIT executable 'run' should exist in ${SH_HOME}"
  exit
else
  echo "                    --> ${SH_HOME}/run exists."
fi

## Delphes ####################################################################################

if [ -z ${DELPHES_HOME} ] ; then
  echo " ######## ERROR ######## "
  echo "In your .cshrc file, please set DELPHES_HOME as the location of the Delphes directory."
  exit
else
  echo "DELPHES_HOME:   ${DELPHES_HOME}"
  if [ ! -d ${DELPHES_HOME} ] ; then 
    echo "${DELPHES_HOME} does not exist. Please create the directory before continuing."
    exit
  fi
fi

DELPHES_VERSION=`grep -e "" ${DELPHES_HOME}/VERSION`
if [ -z ${DELPHES_VERSION} ] ; then
  echo "                  --> Unknown version"
else
  DELPHES_VERSION_PRIMARY=`echo ${DELPHES_VERSION} | sed "s/\([0-9]\).[0-9].[0-9]/\1/g"`
  DELPHES_VERSION_SECONDARY=`echo ${DELPHES_VERSION} | sed "s/[0-9].\([0-9]\).[0-9]/\1/g"`
  DELPHES_VERSION_TERTIARY=`echo ${DELPHES_VERSION} | sed "s/[0-9].[0-9].\([0-9]\)/\1/g"`
  if [ ${DELPHES_VERSION_PRIMARY} != 3 ] ; then
    echo " ######## ERROR ######## "
    echo "** Delphes version 3.X is required for fullscipt.sh to work properly. **"
    exit
  fi
  echo "                    --> v${DELPHES_VERSION}"
fi

if [ ! -e ${DELPHES_HOME}/DelphesSTDHEP ] ; then
  echo " ######## ERROR ######## "
  echo "The Delphes executable 'DelphesSTDHEP' should exist in ${DELPHES_HOME}"
  exit
else
  echo "                    --> ${DELPHES_HOME}/DelphesSTDHEP exists."
fi

## Path and Library System Variables#############################################################

echo ""
echo "LIBPATH:            ${LIBPATH}"
if [[ "${LIBPATH}" == *"${DELPHES_HOME}"* ]] ; then
  echo "                       --> contains '\${DELPHES_HOME}'"
else
  echo " ######## ERROR ######## "
  echo "In your .cshrc file, please add '\${DELPHES_HOME}' to your LIBPATH environment variable."
  exit
fi

echo "LPATH:              ${LPATH}"
if [[ "${LPATH}" == *"${DELPHES_HOME}"* ]] ; then
  echo "                       --> contains '\${DELPHES_HOME}'"
else
  echo " ######## ERROR ######## "
  echo "In your .cshrc file, please add '\${DELPHES_HOME}' to your LPATH environment variable."
  exit
fi

echo "LD_LIBARY_PATH:     ${LD_LIBRARY_PATH}"
if [[ "${LD_LIBRARY_PATH}" == *"${DELPHES_HOME}"* ]] ; then
  echo "                       --> contains '\${DELPHES_HOME}'"
else
  echo " ######## ERROR ######## "
  echo "In your .cshrc file, please add '\${DELPHES_HOME}' to your LD_LIBRARY_PATH environment variable."
  exit
fi

echo "DYLD_LIBRARY_PATH:  ${DYLD_LIBRARY_PATH}"
if [[ "${DYLD_LIBRARY_PATH}" == *"${ROOTSYS}"* ]] ; then
  echo "                       --> contains '\${ROOTSYS}'"
else
  echo " ######## ERROR ######## "
  echo "In your .cshrc file, please add '\${ROOTSYS}' to your DYLD_LIBRARY_PATH environment variable."  exit
  exit
fi

echo "CPLUS_INCLUDE_PATH: ${CPLUS_INCLUDE_PATH}"
if [[ "${CPLUS_INCLUDE_PATH}" == *"${MG_HOME}/MIW_codebase/include"* ]] ; then
  echo "                       --> contains '\${MG_HOME}/MIW_codebase/include'"
else
  echo " ######## ERROR ######## "
  echo "In your .cshrc file, please add '\${MG_HOME}/MIW_codebase/include' to your CPLUS_INCLUDE_PATH environment variable."
  exit
fi
if [[ "${CPLUS_INCLUDE_PATH}" == *"${DELPHES_HOME}"* ]] ; then
  echo "                       --> contains '\${DELPHES_HOME}'"
else
  echo " ######## ERROR ######## "
  echo "In your .cshrc file, please add '\${DELPHES_HOME}' to your CPLUS_INCLUDE_PATH environment variable."
  exit
fi
if [[ "${CPLUS_INCLUDE_PATH}" == *"${DELPHES_HOME}/external"* ]] ; then
  echo "                       --> contains '\${DELPHES_HOME}/external'"
else
  echo " ######## ERROR ######## "
  echo "In your .cshrc file, please add '\${DELPHES_HOME}/external' to your CPLUS_INCLUDE_PATH environment variable."
  exit
fi

echo ""
echo " --> Configuration complete <-- "
echo ""
