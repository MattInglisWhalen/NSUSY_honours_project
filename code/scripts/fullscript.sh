#! /bin/bash

##
## Input: one input file in the style of ${MG_HOME}/MIW_codebase/scripts/input/go_offshell.in
##
## What this script does
## 1) Generates events 
## 2) Analyzes the Delphes output
## 3) Plots the analysis output
##

source ${MG_HOME}/MIW_codebase/scripts/MIW_check_decay_functions.sh
source ${MG_HOME}/MIW_codebase/scripts/MIW_SUSY_constants.sh

## Source the relevant variables
if [ -z "${1}" ] ; then 
  echo "Usage:"
  echo "user> `basename $0` inputfile.in"
  exit
elif [ ! -f ${1} ] ; then
  echo "${1} does not exist"
  echo "Please rerun with a valid input file."
  exit
else
  echo "Sourcing ${1}"
  source ${1}
fi

# Making lists of soft masses that need to be run 
NGRIDPOINTS=0
for (( muMASS=muSTART; muMASS<=muEND; muMASS+=muSTEP )) ; do
for (( M1MASS=M1START; M1MASS<=M1END; M1MASS+=M1STEP )) ; do
for (( M2MASS=M2START; M2MASS<=M2END; M2MASS+=M2STEP )) ; do
for (( M3MASS=M3START; M3MASS<=M3END; M3MASS+=M3STEP )) ; do
for (( Q3MASS=Q3START; Q3MASS<=Q3END; Q3MASS+=Q3STEP )) ; do
for (( U3MASS=U3START; U3MASS<=U3END; U3MASS+=U3STEP )) ; do
  if [ ${DOUBLEOPTWINO} -gt "0" ] ; then
    M2MASS=$((2*M1MASS))
  fi
  if [ ${DOUBLEOPTMU} -gt "0" ] ; then
    muMASS=$((2*M1MASS))
  fi
  muLIST[${NGRIDPOINTS}]="${muMASS}"
  M1LIST[${NGRIDPOINTS}]="${M1MASS}"
  M2LIST[${NGRIDPOINTS}]="${M2MASS}"
  M3LIST[${NGRIDPOINTS}]="${M3MASS}"
  Q3LIST[${NGRIDPOINTS}]="${Q3MASS}"
  U3LIST[${NGRIDPOINTS}]="${U3MASS}"
  NGRIDPOINTS=$((NGRIDPOINTS + 1))
done
done
done
done
done
done

#Then adding the manual ones to the list
SKIPmu=0 ; imu=0
SKIPM1=0 ; iM1=0
SKIPM2=0 ; iM2=0
SKIPM3=0 ; iM3=0
SKIPQ3=0 ; iQ3=0 
SKIPU3=0 ; iU3=0
for muMASS in ${muMAN} ; do
  SKIPmu=0; SKIPM1=0; SKIPM2=0; SKIPM3=0; SKIPQ3=0; SKIPU3=0
for M1MASS in ${M1MAN} ; do
  if [ ${SKIPM1} -lt ${iM1} ] ; then
    SKIPM1=$((SKIPM1+1))
    continue
  fi
for M2MASS in ${M2MAN} ; do
  if [ ${SKIPM2} -lt ${iM2} ] ; then
    SKIPM2=$((SKIPM2+1))
    continue
  fi
for M3MASS in ${M3MAN} ; do
  if [ ${SKIPM3} -lt ${iM3} ] ; then
    SKIPM3=$((SKIPM3+1))
    continue
  fi
for Q3MASS in ${Q3MAN} ; do
  if [ ${SKIPQ3} -lt ${iQ3} ] ; then
    SKIPQ3=$((SKIPQ3+1))
    continue
  fi
for U3MASS in ${U3MAN} ; do
  if [ ${SKIPU3} -lt ${iU3} ] ; then 
    SKIPU3=$((SKIPU3+1))
    continue
  fi
  if [ ${DOUBLEOPTWINO} -gt "0" ] ; then
    M2MASS=$((2*M1MASS))
  fi
  if [ ${DOUBLEOPTMU} -gt "0" ] ; then
    muMASS=$((2*M1MASS))
  fi
  muLIST[${NGRIDPOINTS}]="${muMASS}"
  M1LIST[${NGRIDPOINTS}]="${M1MASS}"
  M2LIST[${NGRIDPOINTS}]="${M2MASS}"
  M3LIST[${NGRIDPOINTS}]="${M3MASS}"
  Q3LIST[${NGRIDPOINTS}]="${Q3MASS}"
  U3LIST[${NGRIDPOINTS}]="${U3MASS}"
  NGRIDPOINTS=$((NGRIDPOINTS + 1))
  break
done ; break
done ; break
done ; break
done ; break
done
  imu=$((imu+1))
  iM1=$((iM1+1))
  iM2=$((iM2+1))
  iM3=$((iM3+1))
  iQ3=$((iQ3+1))
  iU3=$((iU3+1))
done

echo ""
echo "There are ${NGRIDPOINTS} gridpoints in this run"
echo ""

## Make the data storage directories as need be
mkdir -p ${DATA_HOME}/Storage/${STORDIR}/TMODE
mkdir -p ${DATA_HOME}/Storage/${STORDIR}/CMODE
mkdir -p ${DATA_HOME}/Storage/${STORDIR}/BMODE

if [ ${SEWOPT} == 1 ] ; then
  mkdir -p ${DATA_HOME}/Processes/sew_lhe_files
  cp --target-directory=${DATA_HOME}/Processes/sew_lhe_files ${MG_HOME}/MIW_codebase/sew_lhe_files/*.x
fi

if [ ! -d ${MG_HOME}/models/${MODEL} ] ; then
  echo "The model ${MODEL} does not exist in ${MG_HOME}/models"
  echo "Please change the model before rerunning..."
  exit
fi

## Make analysis directories for storage of 
## analyzed .root files and efficiency files
cd ${MG_HOME}
if [ "${ANSYSOPT}" -ge "0" ] ; then
  for MODE in TMODE CMODE BMODE
  do
    for TAG in ${ANSYSTAGS} 
    do
      mkdir -p ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}
    done

    cd ${MG_HOME}/MIW_codebase/include
    NUMATSYS=0
    NUMCMSYS=0
    NUMTOSYS=0
    for TAG in ${ATLTAGS}
    do
      #get out the number of SRs for each ATLAS analysis code
      TEMPNREG=""
      TEMPNREG=`grep -e ".*Int_t nreg_${TAG}=" ${MG_HOME}/MIW_codebase/include/MIW_search_data.h | sed "s/.*Int_t nreg_${TAG}=//g" | sed "s/;//g"`
      if [ -z ${TEMPNREG} ] ; then 
        echo "Analysis ${TAG} does not exist in MIW_search_data.h please remove it from the analysis list"
        exit
      fi
      NREGATL[${NUMATSYS}]=${TEMPNREG}
      NREGTOT[${NUMTOSYS}]=${TEMPNREG}
      if [ -f ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff ] ; then
        echo "The efficiency file ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff already exists..."
        mv ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff_old
      fi
      echo "${TAG} - Number of events = ${NEVENTS}" > ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff
      headerline="${FILEHEADER}"
      for (( j=1;j<=${TEMPNREG};j++ )) 
      do
        headerline="${headerline}, SR${j}"
      done
      echo "${headerline}" >> ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff
      NUMATSYS=`expr ${NUMATSYS} + 1`
      NUMTOSYS=`expr ${NUMTOSYS} + 1`
    done

    for TAG in ${CMSTAGS}
    do
      #get out the number of SRs for each CMS analysis code
      TEMPNREG=`grep -e ".*Int_t nreg_${TAG}=" ${MG_HOME}/MIW_codebase/include/MIW_search_data.h | sed "s/.*Int_t nreg_${TAG}=//g" | sed "s/;//g"`
      NREGCMS[${NUMCMSYS}]=${TEMPNREG}
      NREGTOT[${NUMTOSYS}]=${TEMPNREG}
      if [ -f ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff ] ; then
        echo "The efficiency file ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff already exists..."
        mv ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff_old
      fi
      echo "${TAG} - Number of events = ${NEVENTS}" > ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff
      headerline="${FILEHEADER}"
      for (( j=1;j<=${TEMPNREG};j++ )) 
      do
        headerline="${headerline}, SR${j}"
      done
      echo "${headerline}" >> ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff
      NUMCMSYS=`expr ${NUMCMSYS} + 1`
      NUMTOSYS=`expr ${NUMTOSYS} + 1`
    done
  done
fi

cd ${MG_HOME}
if [ "${GENOPT}" -ge "0" ] ; then

  # Make MadGraph script for ${PROC1DIR} directory creation
  cp ${MG_HOME}/MIW_codebase/templates/mg5_script_template.dat ${MG_HOME}/${PROC1DIR}.dat
  echo "import model ${MODEL}" >> ${MG_HOME}/${PROC1DIR}.dat
  echo "generate ${PROC1GEN1}" >> ${MG_HOME}/${PROC1DIR}.dat
  if [ "${PROC1GEN1}" != "${PROC1GEN2}" ] ; then
    echo "add process ${PROC1GEN2}" >> ${MG_HOME}/${PROC1DIR}.dat
  fi
  echo "output ${DATA_HOME}/Processes/${PROC1DIR}" >> ${MG_HOME}/${PROC1DIR}.dat
  # Run the ${PROC1DIR}.dat script
  ${MG_HOME}/bin/mg5 ${MG_HOME}/${PROC1DIR}.dat
  if [ "${?}" != 0 ] ; then 
    echo ""
    echo " ===> SGE ERROR :  Your chosen core ${HOSTNAME} does not support python 2.6."
    echo "                   Please try running the SGE script again..."
    exit
  fi
  rm ${MG_HOME}/${PROC1DIR}.dat

  # Change the ${PROC1DIR} run_card to nevents/energy required
  sed -i "s/.* = nevents/  ${NEVENTS} = nevents/g" ${DATA_HOME}/Processes/${PROC1DIR}/Cards/run_card.dat
  sed -i "s/.* = ebeam1/  ${CMENERGY1} = ebeam1/g" ${DATA_HOME}/Processes/${PROC1DIR}/Cards/run_card.dat
  sed -i "s/.* = ebeam2/  ${CMENERGY2} = ebeam2/g" ${DATA_HOME}/Processes/${PROC1DIR}/Cards/run_card.dat

  # Make MadGraph script for ${PROC2DIR} directory creation
  cp ${MG_HOME}/MIW_codebase/templates/mg5_script_template.dat ${MG_HOME}/${PROC2DIR}.dat
  echo "import model ${MODEL}" >> ${MG_HOME}/${PROC2DIR}.dat
  echo "generate ${PROC2GEN1}" >> ${MG_HOME}/${PROC2DIR}.dat
  if [ "${PROC2GEN1}" != "${PROC2GEN2}" ] ; then 
    echo "add process ${PROC2GEN2}" >> ${MG_HOME}/${PROC2DIR}.dat
  fi
  echo "output ${DATA_HOME}/Processes/${PROC2DIR}" >> ${MG_HOME}/${PROC2DIR}.dat
  # Run the ${PROC2DIR}.dat script
  ${MG_HOME}/bin/mg5 ${MG_HOME}/${PROC2DIR}.dat
  rm ${MG_HOME}/${PROC2DIR}.dat

  if [ "${SEWOPT}" -eq "1" ] ; then
    NEVENTS=$(( 2 * NEVENTS ))
  fi

  sed -i "s/.* = nevents/  ${NEVENTS} = nevents/g" ${DATA_HOME}/Processes/${PROC2DIR}/Cards/run_card.dat
  sed -i "s/.* = ebeam1/  ${CMENERGY1} = ebeam1/g" ${DATA_HOME}/Processes/${PROC2DIR}/Cards/run_card.dat
  sed -i "s/.* = ebeam2/  ${CMENERGY2} = ebeam2/g" ${DATA_HOME}/Processes/${PROC2DIR}/Cards/run_card.dat

  # Make MadGraph script for ${PROC3DIR} directory creation
  cp ${MG_HOME}/MIW_codebase/templates/mg5_script_template.dat ${MG_HOME}/${PROC3DIR}.dat
  echo "import model ${MODEL}" >> ${MG_HOME}/${PROC3DIR}.dat
  echo "generate ${PROC3GEN1}" >> ${MG_HOME}/${PROC3DIR}.dat
  if [ "${PROC3GEN1}" != "${PROC3GEN2}" ] ; then 
    echo "add process ${PROC3GEN2}" >> ${MG_HOME}/${PROC3DIR}.dat
  fi
  echo "output ${DATA_HOME}/Processes/${PROC3DIR}" >> ${MG_HOME}/${PROC3DIR}.dat
  # Run the ${PROC3DIR}.dat script
  ${MG_HOME}/bin/mg5 ${MG_HOME}/${PROC3DIR}.dat
  rm ${MG_HOME}/${PROC3DIR}.dat

  sed -i "s/.* = nevents/  ${NEVENTS} = nevents/g" ${DATA_HOME}/Processes/${PROC3DIR}/Cards/run_card.dat
  sed -i "s/.* = ebeam1/  ${CMENERGY1} = ebeam1/g" ${DATA_HOME}/Processes/${PROC3DIR}/Cards/run_card.dat
  sed -i "s/.* = ebeam2/  ${CMENERGY2} = ebeam2/g" ${DATA_HOME}/Processes/${PROC3DIR}/Cards/run_card.dat

  ## SUSY-HIT is a public space, so it needs to be replicated
  THIS_SH="${SH_HOME}/../SUSY-HIT_${STORDIR}"
  if [ -d ${THIS_SH} ] && [ ${GENOPT} -ge 0 ] ; then
    ls ${THIS_SH}/.nfs*
    if [ ${?} -eq 0 ] ; then
      echo "This process is already running!"
      echo "End the process then try again..."
      echo ""
      exit
    fi
    rm -r -f ${THIS_SH}
  fi
  cp -r -f ${SH_HOME} ${THIS_SH}
fi

################### RUN SECTION ########################
if [ "${GENOPT}" -ge "0" ] ; then
  ###                                                      ###
  ###         Loop over all necessary gridpoints           ### 
  ###                                                      ###
  for((igridpoint=0;igridpoint<NGRIDPOINTS;igridpoint++))
  do

    muMASS="${muLIST[igridpoint]}"
    M1MASS="${M1LIST[igridpoint]}"
    M2MASS="${M2LIST[igridpoint]}"
    M3MASS="${M3LIST[igridpoint]}"
    Q3MASS="${Q3LIST[igridpoint]}"
    U3MASS="${U3LIST[igridpoint]}"

    MTAG="${muMASS}_${M1MASS}_${M2MASS}_${M3MASS}_${Q3MASS}_${U3MASS}"   #MASSTAG
    echo "Generating events for mu_M1_M2_M3_MQ3_MU3 = ${MTAG}"

    ## Check to make sure this spectrum hasn't already been run
    if [ -d ${DATA_HOME}/Storage/${STORDIR}/TMODE/${MTAG} ] && [ -d ${DATA_HOME}/Storage/${STORDIR}/CMODE/${MTAG} ] && [ -d ${DATA_HOME}/Storage/${STORDIR}/BMODE/${MTAG} ] ; then 
      echo "MTAG ${MTAG} has already been run for all modes, going to next spectrum" 
      continue
    fi

    cd ${MG_HOME}
    ## Copy Les Houches Accord file to SUSY-HIT for width calculations
    cp ${MG_HOME}/MIW_codebase/templates/${MODEL}_susyhit_template.lha ${THIS_SH}/suspect2_lha.in

    cd ${THIS_SH}
    #Changing mu (the higgsino soft mass) 
    sed -i "/Block EXTPAR/,/Block/ s/      ${muTAG}.*/      ${muTAG}     ${muMASS}   # mu(EWSB) (only used if SU_ALGO(6)=0)/g"  ${THIS_SH}/suspect2_lha.in 
    #Changing M_1 (the bino's soft mass)
    sed -i "/Block EXTPAR/,/Block/ s/      ${M1TAG}.*/      ${M1TAG}     ${M1MASS}   # M_1 U(1)_Y Bino mass [~n]  /g" ${THIS_SH}/suspect2_lha.in 
    #Changing M_2 (the wino's soft mass)
    sed -i "/Block EXTPAR/,/Block/ s/      ${M2TAG}.*/      ${M2TAG}     ${M2MASS}   # M_2 SU(2)L Wino mass [~n ~x]  /g" ${THIS_SH}/suspect2_lha.in 
    #Changing M_3 (the gluino's soft mass)
    sed -i "/Block EXTPAR/,/Block/ s/      ${M3TAG}.*/      ${M3TAG}     ${M3MASS}   # M_3 SU(2)c Gluino Mass  /g" ${THIS_SH}/suspect2_lha.in 
    #Changing M_Q3 (the left 3rd gen squark soft mass) 
    sed -i "/Block EXTPAR/,/Block/ s/      ${MQ3TAG}.*/      ${MQ3TAG}     ${Q3MASS}  # M_q3L left 3rd gen squarks/g" ${THIS_SH}/suspect2_lha.in
    #Changing M_U3 (the right stop's soft mass) 
    sed -i "/Block EXTPAR/,/Block/ s/      ${MU3TAG}.*/      ${MU3TAG}     ${U3MASS}   # M_U3  right stop /g" ${THIS_SH}/suspect2_lha.in
    ${THIS_SH}/run  # you must make sure that susyhit is set to run with SuSpect
    rm -f ${THIS_SH}/sed*

    ## This extraction method relies on the mass block occurring before the decay block
    #get out the neutralino mass
    grep -e "   ${N1LABEL}" ${THIS_SH}/susyhit_slha.out | sed "s/   ${N1LABEL}    //g" | sed "s/  #.*//g" > ${THIS_SH}/${N1LABEL}.mass
    exec<${THIS_SH}/${N1LABEL}.mass
    read line
    TEMPN1=${line}
    #get out the chargino mass
    grep -e "   ${X1LABEL}" ${THIS_SH}/susyhit_slha.out | sed "s/   ${X1LABEL}    //g" | sed "s/  #.*//g" > ${THIS_SH}/${X1LABEL}.mass
    exec<${THIS_SH}/${X1LABEL}.mass
    read line
    TEMPX1=${line}
    #get out the gluino mass
    grep -e "   ${GOLABEL}" ${THIS_SH}/susyhit_slha.out | sed "s/   ${GOLABEL}    //g" | sed "s/  #.*//g" > ${THIS_SH}/${GOLABEL}.mass
    exec<${THIS_SH}/${GOLABEL}.mass
    read line
    TEMPGO=${line}
    #get out the stop 1 mass
    grep -e "   ${T1LABEL}" ${THIS_SH}/susyhit_slha.out | sed "s/   ${T1LABEL}    //g" | sed "s/  #.*//g" > ${THIS_SH}/${T1LABEL}.mass
    exec<${THIS_SH}/${T1LABEL}.mass
    read line
    TEMPT1=${line}
    rm -f ${THIS_SH}/*.mass

    echo "The hard masses are    ${TEMPN1} ${TEMPX1} ${TEMPGO} ${TEMPT1}"
    echo "As integers, these are `${MG_HOME}/MIW_codebase/bin/get_int_mass_from_SH.x ${TEMPN1}`               \
`${MG_HOME}/MIW_codebase/bin/get_int_mass_from_SH.x ${TEMPX1}`                `${MG_HOME}/MIW_codebase/bin/get_int_mass_from_SH.x ${TEMPGO}`\
                `${MG_HOME}/MIW_codebase/bin/get_int_mass_from_SH.x ${TEMPT1}` "

    allowed_by_opt ${DECOPT1} ${TEMPN1} ${TEMPX1} ${TEMPGO} ${TEMPT1}
    ALLOWED1=$?
    allowed_by_opt ${DECOPTT} ${TEMPN1} ${TEMPX1} ${TEMPGO} ${TEMPT1}
    ALLOWEDT=$?
    allowed_by_opt ${DECOPTB} ${TEMPN1} ${TEMPX1} ${TEMPGO} ${TEMPT1}
    ALLOWEDB=$?

    echo "Allowed flags: ${ALLOWED1} ${ALLOWEDT} ${ALLOWEDB}"
    
    #in this process the gluino needs to be able to decay to a stop/top pair
    if [ ${ALLOWED1} -le 0 ] ; then 
      echo "The decay process ${DECOPT1} (see MIW_script_functions.h) is kinematically forbidden, trying the next spectrum..."
      continue
    fi

    ## Check to make sure this spectrum hasn't already been run
    if [ ${ALLOWEDT} -gt "0" ] && [ ! -d ${DATA_HOME}/Storage/${STORDIR}/TMODE/${MTAG} ] ; then 
      echo "T-mode allowed but directory TMODE/${MTAG} doesn't exist, running it" 
    else
      if [ ${ALLOWEDB} -gt "0" ] && [ ! -d ${DATA_HOME}/Storage/${STORDIR}/BMODE/${MTAG} ] ; then 
        echo "B-mode allowed but directory BMODE/${MTAG} doesn't exist, running it" 
      else
        if [ ${ALLOWEDT} -gt "0" ] && [ ${ALLOWEDB} -gt "0" ] && [ ! -d ${DATA_HOME}/Storage/${STORDIR}/CMODE/${MTAG} ] ; then 
          echo "C-mode allowed but directory CMODE/${MTAG} doesn't exist, running it" 
        else
          continue  
        fi
      fi
    fi

    cp ${THIS_SH}/susyhit_slha.out ${DATA_HOME}/Processes/${PROC1DIR}/Cards/param_card.dat
    cp ${THIS_SH}/susyhit_slha.out ${DATA_HOME}/Processes/${PROC2DIR}/Cards/param_card.dat
    cp ${THIS_SH}/susyhit_slha.out ${DATA_HOME}/Processes/${PROC3DIR}/Cards/param_card.dat

    date
    for DIR in ${PROC1DIR} ${PROC2DIR} ${PROC3DIR}
    do
      if [ "${DIR}" == "${PROC2DIR}" ] && [ -d ${DATA_HOME}/Storage/${STORDIR}/CMODE/${MTAG} ] && [ -d ${DATA_HOME}/Storage/${STORDIR}/TMODE/${MTAG} ] ; then continue ; fi # already have tmode events
      if [ "${DIR}" == "${PROC3DIR}" ] && [ -d ${DATA_HOME}/Storage/${STORDIR}/CMODE/${MTAG} ] && [ -d ${DATA_HOME}/Storage/${STORDIR}/BMODE/${MTAG} ] ; then continue ; fi # already have bmode events
      if [ "${SEWOPT}" == 0 ] && [ "${DIR}" == "${PROC1DIR}" ] && [ -d ${DATA_HOME}/Storage/${STORDIR}/CMODE/${MTAG} ] ; then continue ; fi # already have cmode events, but sewing needs the production events
      if [ "${DIR}" == "${PROC2DIR}" ] && [ ${ALLOWEDT} -le 0 ] ; then continue ; fi  # Tmode is kinematically forbidden
      if [ "${DIR}" == "${PROC3DIR}" ] && [ ${ALLOWEDB} -le 0 ] ; then continue ; fi  # Bmode is kinematically forbidden
      if [ "${SEWOPT}" == 0 ] && [ "${DIR}" == "${PROC1DIR}" ] ; then if [ ${ALLOWEDT} -le 0 ] || [ ${ALLOWEDB} -le 0 ] ; then continue ; fi ; fi # Cmode is kinematically forbidden
      cd ${MG_HOME}
      rm -r -f ${DATA_HOME}/Processes/${DIR}/Events/*
      cp ${MG_HOME}/MIW_codebase/templates/mg5_script_template.dat ${MG_HOME}/${DIR}.dat
      echo "launch ${DATA_HOME}/Processes/${DIR} -i" >> ${MG_HOME}/${DIR}.dat
      if [ "${SEWOPT}" -eq "1" ] ; then  ## we need to sew these together first
        echo "generate_events --laststep=parton -f" >> ${MG_HOME}/${DIR}.dat
      else                               ## go right ahead and shower
        echo "generate_events --laststep=pythia -f" >> ${MG_HOME}/${DIR}.dat
      fi
      # Everything is set up, let's run the production process
      ${MG_HOME}/bin/mg5 ${MG_HOME}/${DIR}.dat
    done

    date
    ##this is where we sew the production events together with the decay events
    ## We are in a public space so everything needs to be tagged as _${STORDIR}
    if [ "${SEWOPT}" -eq "1" ] ; then 
      cd ${DATA_HOME}/Processes/sew_lhe_files
      mv -f ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/unweighted_events.lhe.gz ${DATA_HOME}/Processes/sew_lhe_files/unweighted_events_prod_${STORDIR}.lhe.gz 
      if [ ${ALLOWEDT} -gt 0 ] && [ -f ${DATA_HOME}/Processes/${PROC2DIR}/Events/run_01/unweighted_events.lhe.gz ] ; then 
        mv ${DATA_HOME}/Processes/${PROC2DIR}/Events/run_01/unweighted_events.lhe.gz ${DATA_HOME}/Processes/sew_lhe_files/unweighted_events_tmode_${STORDIR}.lhe.gz
      fi
      if [ ${ALLOWEDB} -gt 0 ] && [ -f ${DATA_HOME}/Processes/${PROC3DIR}/Events/run_01/unweighted_events.lhe.gz ] ; then 
        mv ${DATA_HOME}/Processes/${PROC3DIR}/Events/run_01/unweighted_events.lhe.gz ${DATA_HOME}/Processes/sew_lhe_files/unweighted_events_bmode_${STORDIR}.lhe.gz 
      fi
      rm -f *_${STORDIR}.lhe
      gzip -d *_${STORDIR}.lhe.gz
      if [ "${?}" != 0 ] ; then echo "I don't know what went wrong..." ; fi
      rm -f *_${STORDIR}.in
      ##TMODE
      if [ ${ALLOWEDT} -gt 0 ] ; then 
        if [ -f unweighted_events_tmode_${STORDIR}.lhe ] ; then
          echo "-------- Sewing the T-mode files together ---------"
          echo "${PRODPART}" > sew_tmode_${STORDIR}.in
          echo "unweighted_events_prod_${STORDIR}.lhe" >> sew_tmode_${STORDIR}.in
          echo "unweighted_events_tmode_${STORDIR}.lhe" >> sew_tmode_${STORDIR}.in
          echo "tmode_${STORDIR}.out" >> sew_tmode_${STORDIR}.in
          ./sew_full.x < sew_tmode_${STORDIR}.in
          if [ "`echo $?`" -ne 0 ] ; then echo "./sew_full.x failure, trying next spectrum..." ; continue ; fi
          nt=`wc -l tmode_${STORDIR}.out | sed "s/ tmode_${STORDIR}.out//g"`
          head -$((nt-7)) tmode_${STORDIR}.out > sewed_unweighted_events_tmode_${STORDIR}.lhe
          gzip sewed_unweighted_events_tmode_${STORDIR}.lhe
          mv sewed_unweighted_events_tmode_${STORDIR}.lhe.gz ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/unweighted_events_TMODE.lhe.gz
        fi
        ##CMODE
        if [ ${ALLOWEDB} -gt 0 ] ; then 
          if [ -f unweighted_events_bmode_${STORDIR}.lhe ] ; then
            echo ""
            echo "-------- Sewing the C-mode files together : Step 1 ---------"
            echo "${PRODPART}" > sew_cmode1_${STORDIR}.in
            echo "unweighted_events_prod_${STORDIR}.lhe" >> sew_cmode1_${STORDIR}.in
            echo "unweighted_events_tmode_${STORDIR}.lhe" >> sew_cmode1_${STORDIR}.in
            echo "cmode1_${STORDIR}.out" >> sew_cmode1_${STORDIR}.in
            ./sew_half.x < sew_cmode1_${STORDIR}.in
            if [ "`echo $?`" -ne 0 ] ; then echo "./sew_half.x failure, trying next spectrum..." ; continue ; fi
            nc1=`wc -l cmode1_${STORDIR}.out | sed "s/ cmode1_${STORDIR}.out//g"`
            head -$((nc1-7)) cmode1_${STORDIR}.out > cmode2_${STORDIR}.out
            echo "-------- Sewing the C-mode files together : Step 2 ---------"
            echo "${GOLABEL}" > sew_cmode3_${STORDIR}.in
            echo "cmode2_${STORDIR}.out" >> sew_cmode3_${STORDIR}.in
            echo "unweighted_events_bmode_${STORDIR}.lhe" >> sew_cmode3_${STORDIR}.in
            echo "cmode3_${STORDIR}.out" >> sew_cmode3_${STORDIR}.in
            ./sew_half.x < sew_cmode3_${STORDIR}.in
            if [ "`echo $?`" -ne 0 ] ; then echo "./sew_half.x failure, trying next spectrum..." ; continue ; fi
            nc3=`wc -l cmode3_${STORDIR}.out | sed "s/ cmode3_${STORDIR}.out//g"`
            head -$((nc3-7)) cmode3_${STORDIR}.out > sewed_unweighted_events_cmode_${STORDIR}.lhe
            gzip sewed_unweighted_events_cmode_${STORDIR}.lhe
            mv sewed_unweighted_events_cmode_${STORDIR}.lhe.gz ../${PROC1DIR}/Events/run_01/unweighted_events_CMODE.lhe.gz
          fi
        fi
      fi
      ##BMODE
      if [ ${ALLOWEDB} -gt 0 ] ; then
        if [ -f unweighted_events_bmode_${STORDIR}.lhe  ] ; then
          echo ""
          echo "-------- Sewing the B-mode files together ---------"
          echo "${PRODPART}" > sew_bmode_${STORDIR}.in
          echo "unweighted_events_prod_${STORDIR}.lhe" >> sew_bmode_${STORDIR}.in
          echo "unweighted_events_bmode_${STORDIR}.lhe" >> sew_bmode_${STORDIR}.in
          echo "bmode_${STORDIR}.out" >> sew_bmode_${STORDIR}.in
          ./sew_full.x < sew_bmode_${STORDIR}.in
          if [ "`echo $?`" -ne 0 ] ; then echo "./sew_full.x failure, trying next spectrum..." ; continue ; fi
          nb=`wc -l bmode_${STORDIR}.out | sed "s/ bmode_${STORDIR}.out//g"`
          head -$((nb-7)) bmode_${STORDIR}.out > sewed_unweighted_events_bmode_${STORDIR}.lhe
          gzip sewed_unweighted_events_bmode_${STORDIR}.lhe
          mv sewed_unweighted_events_bmode_${STORDIR}.lhe.gz ../${PROC1DIR}/Events/run_01/unweighted_events_BMODE.lhe.gz
        fi
      fi
      rm -f *_${STORDIR}.lhe *_${STORDIR}.out

      date
      for MODE in TMODE CMODE BMODE
      do
        cd ${MG_HOME}
        if [ -d ${DATA_HOME}/Storage/${STORDIR}/${MODE}/${MTAG} ] ; then continue ; fi
        if [ ${ALLOWEDT} -le 0 ] ; then if [ "${MODE}" == "TMODE" -o "${MODE}" == "CMODE" ] ; then continue ; fi ; fi
        if [ ${ALLOWEDB} -le 0 ] ; then if [ "${MODE}" == "BMODE" -o "${MODE}" == "CMODE" ] ; then continue ; fi ; fi
        mv ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/unweighted_events_${MODE}.lhe.gz ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/unweighted_events.lhe.gz 
        #Make the MadGraph script to run the sewn together process
        cp ${MG_HOME}/MIW_codebase/templates/mg5_script_template.dat ${MG_HOME}/${STORDIR}_${MODE}.dat
        echo "launch ${DATA_HOME}/Processes/${PROC1DIR} -i" >> ${MG_HOME}/${STORDIR}_${MODE}.dat
        echo "pythia run_01 --tag=${MODE} --laststep=pythia -f" >> ${MG_HOME}/${STORDIR}_${MODE}.dat
        # Everything is set up, let's run MadGraph5
        ${MG_HOME}/bin/mg5 ${MG_HOME}/${STORDIR}_${MODE}.dat
        rm -f ${MG_HOME}/${STORDIR}_${MODE}.dat
        if [ ! -f ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/${MODE}_pythia_events.hep.gz ] ; then
          echo "Failed to create pythia files for MTAG ${MTAG} with mode ${MODE}"
          continue
        fi
        mkdir -p ${DATA_HOME}/Storage/${STORDIR}/${MODE}/${MTAG}
        cp ${DATA_HOME}/Processes/${PROC1DIR}/Cards/param_card.dat ${DATA_HOME}/Storage/${STORDIR}/${MODE}/${MTAG}/param_card.dat #For future reference
      done
    else # SEWOPT = 0
      mkdir -p ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01
      if [ ${ALLOWEDT} -gt 0 ] ; then
        mkdir -p ${DATA_HOME}/Storage/${STORDIR}/TMODE/${MTAG}
        cp ${DATA_HOME}/Processes/${PROC1DIR}/Cards/param_card.dat ${DATA_HOME}/Storage/${STORDIR}/TMODE/${MTAG}/param_card.dat #For future reference
        mv ${DATA_HOME}/Processes/${PROC2DIR}/Events/run_01/tag_1_pythia_events.hep.gz ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/TMODE_pythia_events.hep.gz
        if [ ${ALLOWEDB} -gt 0 ] ; then
          mkdir -p ${DATA_HOME}/Storage/${STORDIR}/CMODE/${MTAG}
          cp ${DATA_HOME}/Processes/${PROC1DIR}/Cards/param_card.dat ${DATA_HOME}/Storage/${STORDIR}/CMODE/${MTAG}/param_card.dat #For future reference
          if [ -f ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/tag_1_pythia_events.hep.gz ] ; then
            mv ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/tag_1_pythia_events.hep.gz ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/CMODE_pythia_events.hep.gz
          fi
        fi
      fi
      if [ ${ALLOWEDB} -gt 0 ] ; then
        mkdir -p ${DATA_HOME}/Storage/${STORDIR}/BMODE/${MTAG}
        cp ${DATA_HOME}/Processes/${PROC1DIR}/Cards/param_card.dat ${DATA_HOME}/Storage/${STORDIR}/BMODE/${MTAG}/param_card.dat #For future reference
        if [ -f ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/tag_1_pythia_events.hep.gz ] ; then
          mv ${DATA_HOME}/Processes/${PROC3DIR}/Events/run_01/tag_1_pythia_events.hep.gz ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/BMODE_pythia_events.hep.gz
        fi
      fi
    fi # if SEWOPT=1

    for MODE in TMODE CMODE BMODE 
    do
      if [ ! -f ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/${MODE}_pythia_events.hep.gz ] ; then 
        echo "${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/${MODE}_pythia_events.hep.gz doesn't exist"
        echo "Maybe this spectrum/mode was already run. Trying next mode"
        continue
      fi
      gzip -d ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/${MODE}_pythia_events.hep.gz
      if [ "${GENOPT}" -ne "2" ] ; then
        if [ -f ${DATA_HOME}/Storage/${STORDIR}/${MODE}/${MTAG}/tag_1_delphes_events_ATLAS.root ] ; then
          echo "${DATA_HOME}/Storage/${STORDIR}/${MODE}/${MTAG}/tag_1_delphes_events_ATLAS.root already exists?"
        else
          ${DELPHES_HOME}/DelphesSTDHEP ${MG_HOME}/Template/Cards/delphes_card_ATLAS.dat ${DATA_HOME}/Storage/${STORDIR}/${MODE}/${MTAG}/tag_1_delphes_events_ATLAS.root ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/${MODE}_pythia_events.hep
          cd ${DATA_HOME}/Storage/${STORDIR}/${MODE}/${MTAG}
          cp --target-directory=./ ${MG_HOME}/MIW_codebase/bin/run_root_script_with_delphes.C ${MG_HOME}/MIW_codebase/bin/run_trim_delphes.C
          root -l -q -b run_root_script_with_delphes.C++\(\"run_trim_delphes\",\"tag_1_delphes_events_ATLAS.root\"\)
          if [ "$?" -eq 0 ] ; then
            mv trimmed_delphes.root tag_1_delphes_events_ATLAS.root
          else 
            echo "ERROR while trimming delphes.root"
          fi
          rm -f run_*
        fi
      fi
      if [ "${GENOPT}" -ne "1" ] ; then
        if [ -f ${DATA_HOME}/Storage/${STORDIR}/${MODE}/${MTAG}/tag_1_delphes_events_CMS.root ] ; then
          echo "${DATA_HOME}/Storage/${STORDIR}/${MODE}/${MTAG}/tag_1_delphes_events_CMS.root already exists?"
        else
          ${DELPHES_HOME}/DelphesSTDHEP ${MG_HOME}/Template/Cards/delphes_card_CMS.dat ${DATA_HOME}/Storage/${STORDIR}/${MODE}/${MTAG}/tag_1_delphes_events_CMS.root ${DATA_HOME}/Processes/${PROC1DIR}/Events/run_01/${MODE}_pythia_events.hep
          cd ${DATA_HOME}/Storage/${STORDIR}/${MODE}/${MTAG}
          cp --target-directory=./ ${MG_HOME}/MIW_codebase/bin/run_root_script_with_delphes.C ${MG_HOME}/MIW_codebase/bin/run_trim_delphes.C
          root -l -q -b run_root_script_with_delphes.C++\(\"run_trim_delphes\",\"tag_1_delphes_events_CMS.root\"\)
          if [ "$?" -eq 0 ] ; then
            mv trimmed_delphes.root tag_1_delphes_events_CMS.root
          else 
            echo "ERROR while trimming delphes.root"
          fi
          rm -f run_*
        fi
      fi
    done

    rm -r -f ${DATA_HOME}/Processes/${PROC1DIR}/Events/*
    rm -r -f ${DATA_HOME}/Processes/${PROC2DIR}/Events/*
    rm -r -f ${DATA_HOME}/Processes/${PROC3DIR}/Events/*

  done # for loop over gridpoints

  rm -f ${MH_HOME}/${PROC1DIR}.dat ${MH_HOME}/${PROC2DIR}.dat ${MH_HOME}/${PROC3DIR}.dat

fi  ## genopt -ge 0
    
################ ANALYSIS SECTION ###################
if [ "${ANSYSOPT}" -ge "0" ] ; then
  for((igridpoint=0;igridpoint<NGRIDPOINTS;igridpoint++))
  do

    muMASS="${muLIST[igridpoint]}"
    M1MASS="${M1LIST[igridpoint]}"
    M2MASS="${M2LIST[igridpoint]}"
    M3MASS="${M3LIST[igridpoint]}"
    Q3MASS="${Q3LIST[igridpoint]}"
    U3MASS="${U3LIST[igridpoint]}"
    if [ -z ${muMASS} ] ; then echo "muMASS - something went horribly wrong" ; igridpoint=$((igridpoint-1)) ; continue ; fi
    if [ -z ${M1MASS} ] ; then echo "M1MASS - something went horribly wrong" ; igridpoint=$((igridpoint-1)) ; continue ; fi
    if [ -z ${M2MASS} ] ; then echo "M2MASS - something went horribly wrong" ; igridpoint=$((igridpoint-1)) ; continue ; fi
    if [ -z ${M3MASS} ] ; then echo "M3MASS - something went horribly wrong" ; igridpoint=$((igridpoint-1)) ; continue ; fi
    if [ -z ${Q3MASS} ] ; then echo "Q3MASS - something went horribly wrong" ; igridpoint=$((igridpoint-1)) ; continue ; fi
    if [ -z ${U3MASS} ] ; then echo "U3MASS - something went horribly wrong" ; igridpoint=$((igridpoint-1)) ; continue ; fi
    MTAG="${muMASS}_${M1MASS}_${M2MASS}_${M3MASS}_${Q3MASS}_${U3MASS}"   #MASSTAG
    echo "Analyzing events for mu_M1_M2_M3_MQ3_MU3 = ${MTAG}"

    for MODE in TMODE CMODE BMODE
    do
      cd ${DATA_HOME}/Storage/${STORDIR}/${MODE}
      if [ -d ${MTAG} ] ; then
        cd ${MTAG}
      else
        echo "${MTAG} spectrum has not been run for ${MODE}, trying next mode/spectrum..."
        continue
      fi

      ## Placing particle masses into variables for later
      ## Masses are kept in the param_card.dat in the event directories
      ## This extraction method relies on the mass block occurring before the decay block
      #get out the neutralino mass
      grep -e "   ${N1LABEL}" param_card.dat | sed "s/   ${N1LABEL}    //g" | sed "s/  #.*//g" > ${N1LABEL}.mass
      exec<${N1LABEL}.mass
      read line
      TEMPN1=${line}
      #get out the chargino mass
      grep -e "   ${X1LABEL}" param_card.dat | sed "s/   ${X1LABEL}    //g" | sed "s/  #.*//g" > ${X1LABEL}.mass
      exec<${X1LABEL}.mass
      read line
      TEMPX1=${line}
      #get out the gluino mass
      grep -e "   ${GOLABEL}" param_card.dat | sed "s/   ${GOLABEL}    //g" | sed "s/  #.*//g" > ${GOLABEL}.mass
      exec<${GOLABEL}.mass
      read line
      TEMPGO=${line}
      #get out the sbottom 1 mass
      grep -e "   ${B1LABEL}" param_card.dat | sed "s/   ${B1LABEL}    //g" | sed "s/  #.*//g" > ${B1LABEL}.mass
      exec<${B1LABEL}.mass
      read line
      TEMPB1=${line}
      #get out the stop 1 mass
      grep -e "   ${T1LABEL}" param_card.dat | sed "s/   ${T1LABEL}    //g" | sed "s/  #.*//g" > ${T1LABEL}.mass
      exec<${T1LABEL}.mass
      read line
      TEMPT1=${line}
      LINELABEL="${muMASS}, ${M1MASS}, ${M2MASS}, ${M3MASS}, ${Q3MASS}, ${U3MASS},${TEMPN1},${TEMPX1},${TEMPGO},${TEMPB1},${TEMPT1}"
      rm -f *.mass sed*

      cp --target-directory=./ ${MG_HOME}/MIW_codebase/bin/run_root_script_with_delphes.C ${MG_HOME}/MIW_codebase/bin/run_root_analysis.C
      # -- ATLAS analysis -- #
      if [ "${ANSYSOPT}" -ne "2" ] ; then
        if [ -f tag_1_delphes_events_ATLAS.root ] ; then
          for TAG in ${ATLTAGS}
          do
            if [ ! -f "${TAG}_eff.txt" ] ; then
              root -l -q -b run_root_script_with_delphes.C++\(\"run_root_analysis\",\"${TAG}\"\)  # run it
              mv -f ${TAG}_analysis.root ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${MTAG}.root
            else
              echo "${TAG}_eff.txt already exists so we won't redo the analysis"
            fi
            echo "Placing ${LINELABEL}, `grep -e "" ${TAG}_eff.txt` into efficiency file."
            echo "${LINELABEL}, `grep -e "" ${TAG}_eff.txt`" >> ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff
          done
        else
          echo "Delphes files were not generated..."
          it=0
          for TAG in ${ATLTAGS}
          do
            fillline=""
            for (( j=1;j<=${NREGATL[it]};j++ )) 
            do
              fillline="0, ${fillline}"
            done
            rm -f ${TAG}_eff.temp
            echo "${fillline}" > ${TAG}_eff.temp
            it=$((it + 1))
            echo "Placing ${LINELABEL}, `grep -e "" ${TAG}_eff.temp` into efficiency file."
            echo "${LINELABEL}, `grep -e "" ${TAG}_eff.temp`" >> ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff
          done
        fi
      fi  ##ansysopt -ne 2

      # -- CMS ANALYSIS -- #
      if [ "${ANSYSOPT}" -ne "1" ] ; then
        if [ -f tag_1_delphes_events_CMS.root ] ; then
          for TAG in ${CMSTAGS}
          do
            if [ ! -f "${TAG}_eff.txt" ] ; then
              root -l -q -b run_root_script_with_delphes.C++\(\"run_root_analysis\",\"${TAG}\"\)  # run it
              mv -f ${TAG}_analysis.root ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${MTAG}.root
            else
              echo "${TAG}_eff.txt already exists so we won't redo the analysis"
            fi
            echo "Placing ${LINELABEL}, `grep -e "" ${TAG}_eff.txt` into efficiency file."
            echo "${LINELABEL}, `grep -e "" ${TAG}_eff.txt`" >> ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff
          done
        else
          echo "Delphes files were not generated..."
          for TAG in ${CMSTAGS}
          do
            fillline=""
            for (( j=1;j<=NREGCMS[it];j++ )) 
            do
              fillline="0, ${fillline}"
            done
            rm -f ${TAG}_eff.temp
            echo "${fillline}" > ${TAG}_eff.temp
            it=$((it + 1))
            echo "Placing ${LINELABEL}, `grep -e "" ${TAG}_eff.temp` into temporary efficiency file."
            echo "${LINELABEL}, `grep -e "" ${TAG}_eff.temp`" >> ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/${TAG}.eff
          done
        fi
      fi  ##ansysopt -ne 1
      rm -f run_*
    done # modes

    #Done!
    echo "Done analyzing ${MTAG} !!"

  done
fi  ##ansysopt -ge 0

date
if [ "${PLOTOPT}" -ge "0" ] ; then
  ## Splitting the efficiency file into mu_M1_M2 pieces for plotting purposes
  cd ${MG_HOME}
  CURRTAG="SCRATCH"
  for MODE in TMODE CMODE BMODE
  do
    for TAG in ${ANSYSTAGS} 
    do
      cd ${MG_HOME}/analysis/${STORDIR}_${MODE}/${TAG}/
      rm -f *_${TAG}.eff
      exec<${TAG}.eff
      read TITLELINE
      read COLHEADER
      while read line
      do
        THISTAG=` echo ${line} | sed "s/\([0-9]*, [0-9]*, [0-9]*\).*/\1/g" | sed "s/, /_/g" `
        if [ ! -f "${THISTAG}_${TAG}.eff" ] ; then
          echo "${TITLELINE}" > ${THISTAG}_${TAG}.eff
          echo "${COLHEADER}" >> ${THISTAG}_${TAG}.eff
        fi
        echo "${line}" >> ${THISTAG}_${TAG}.eff
      done
    done
  done

  mkdir -p ${MG_HOME}/analysis/${STORDIR}

  isys=0
  for TAG in ${ANSYSTAGS}
  do
    cd ${MG_HOME}/analysis/${STORDIR}_TMODE/${TAG}
    for EFFFILE in ` ls *${TAG}.eff* `
    do
      ADD=` echo ${EFFFILE} | sed "s/_${TAG}.eff//g" | sed "s/${TAG}.eff//g"  | sed "s/_old//g" `
      UPDTAGS[${isys}]="${UPDTAGS[isys]} ${ADD}"
    done
    isys=$((isys + 1))
  done

  isys=0
  cd ${MG_HOME}/analysis
  cp ${MG_HOME}/MIW_codebase/plotting/plot_excl_br_curves_${STORDIR}.C ./
  for TAG in ${ANSYSTAGS}
  do
    for UPD in ${UPDTAGS[isys]}
    do
      echo "${UPD}"
      root -l -q -b plot_excl_br_curves_${STORDIR}.C++\(\"${STORDIR}\",\"${TAG}\",\"${UPD}\"\)
      mv --target-directory=${STORDIR} *${TAG}_${STORDIR}*.pdf
      rm ${TAG}_exclusion.root
      if [ "${STORDIR}" == "GO_OFFSHELL" ] ; then break ; fi
      if [ "${STORDIR}" == "SCAN_T1PAIR" ] ; then break ; fi
    done
    isys=$((isys + 1))
  done

  rm plot_excl_br_curves_${STORDIR}*

fi ##PLOTOPT ge 0

date

