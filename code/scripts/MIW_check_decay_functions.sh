#!/bin/bash

TOPMASS=175
TOPMASS2=$((2 * TOPMASS))
BOTMASS=5
BOTMASS2=$((2 * BOTMASS))

# for gluino pair production with go>ttxn1 through an offshell stop
function go_tmode_allowed { ##arg1==gluino mass, arg2==neut mass
  GOMASS=`get_int_mass_from_SH.x ${1}`
  N1MASS=`get_int_mass_from_SH.x ${2}`
  if [ $((GOMASS)) -gt $((N1MASS+TOPMASS2)) ] ; then
    return 1
  else
    return 0
  fi
}

# for gluino pair production with go>bbxn1 through an offshell sbottom
function go_bmode_allowed { ##arg1==gluino mass, arg2==neut mass
  GOMASS=`get_int_mass_from_SH.x ${1}`
  N1MASS=`get_int_mass_from_SH.x ${2}`
  if [ $((GOMASS)) -gt $((N1MASS+BOTMASS2)) ] ; then
    return 1
  else
    return 0
  fi
}

#for gluinos decaying to a stop-top pair
function go_t1t_allowed { ## arg1==gluino mass, arg2==stop mass
  GOMASS=`get_int_mass_from_SH.x ${1}`
  T1MASS=`get_int_mass_from_SH.x ${2}`
  if [ $((GOMASS)) -gt $((T1MASS+TOPMASS)) ]; then
    return 1
  else
    return 0
  fi
}

#for on shell stops decaying to a top neutralino
function t1_tmode_allowed { ## arg1== stop mass, arg2==neut mass
  T1MASS=`get_int_mass_from_SH.x ${1}`
  N1MASS=`get_int_mass_from_SH.x ${2}`
  if [ $((T1MASS)) -gt $((N1MASS+TOPMASS)) ]; then
    return 1
  else
    return 0
  fi
}

#for on shell stops decaying to a bottom chargino
function t1_bmode_allowed { ## arg1==stop mass, arg2==chrg mass
  T1MASS=`get_int_mass_from_SH.x ${1}`
  X1MASS=`get_int_mass_from_SH.x ${2}`
  if [ $((T1MASS)) -gt $((X1MASS+BOTMASS)) ]; then
    return 1
  else
    return 0
  fi
}

## arg1==option arg2==neut1 arg3==chrg1 arg4==gluino arg4==stop1
function allowed_by_opt { 
  OPT=${1}
  fN1=${2}
  fX1=${3}
  fGO=${4}
  fT1=${5}
  if [ "${OPT}" -eq "0" ] ; then
    return 1
  elif [ "${OPT}" -eq "1" ] ; then
    go_tmode_allowed ${fGO} ${fN1}
    return $?
  elif [ "${OPT}" -eq "2" ] ; then 
    go_bmode_allowed ${fGO} ${fN1}
    return $?
  elif [ "${OPT}" -eq "3" ] ; then 
    go_t1t_allowed ${fGO} ${fT1}
    return $?
  elif [ "${OPT}" -eq "4" ] ; then 
    t1_tmode_allowed ${fT1} ${fN1}
    return $?
  elif [ "${OPT}" -eq "5" ] ; then 
    t1_bmode_allowed ${fT1} ${fX1}
    return $?
  else
    echo "allowed_by_opt: option ${OPT} not allowed"
    return 0
  fi
}