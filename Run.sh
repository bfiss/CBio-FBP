#!/bin/bash
# This script should be at rosetta/targets/

E_BADARGS=65

function echo_stage_bool_opt {
  args=("$@")
  for i in 1 2 3 4 5
  do
    if [ ${args[i-1]} -eq "0" ]
    then
      echo "-compbio_app:skip_stage${i} false"
    else
      echo "-compbio_app:skip_stage${i} true"
    fi 
  done
}

function echo_stage_num_opt {
  args=("$@")
  for i in 1 2 3 4 5
  do
    echo "-compbio_app:num_cycles_stage${i} ${args[i-1]}"
  done
}

function echo_frags_opt {
  if [ $1 -eq 0 ]
  then
    echo "-compbio_app:large_frags false"
  else
    echo "-compbio_app:large_frags true"
  fi
}

function echo_first_run_opt {
  if [ $1 -eq 0 ]
  then
    echo "-compbio_app:first_run false"
    echo "-compbio_app:last_structure model_c.pdb"
  else
    echo "-compbio_app:first_run true"
  fi
}

#
# Arguments:
# 
# 1 = temperature
# 2-6 = if stage# should be skipped
# 7-11 = number of cycles of the stage#
# 12 = if this is the first run
# 13 = if large fragments should be used
# 14 = current number of aminoacids in the protein
# 15 = name of the protein
# 16 = number of outputs to be generated
#
# Result: outputs the contents of the desired flag file
#

function echo_options {
  echo "-in:file:fasta inputs/${15}.fasta1"
  echo "-in:file:frag3 inputs/fragment_files/aat`printf "%03d" ${14}`_03_05.200_v1_3"
  echo "-in:file:frag9 inputs/fragment_files/aat`printf "%03d" ${14}`_09_05.200_v1_3"
  echo "-out:nstruct ${16}"
  echo "-out:pdb true"
  echo "-mute core.io.database"
  echo "-compbio_app:temperature ${1}"
  echo_stage_bool_opt $2 $3 $4 $5 $6
  echo_stage_num_opt $7 $8 $9 ${10} ${11}
  echo_first_run_opt ${12}
  echo_frags_opt ${13}
  echo ""
}

#
# Arguments:
#
# 1 the FASTA sequence
# 2 the current protein size
# 3 number of the trial
# 4-17 echo_options arguments
# 18 name of the protein
#
# Does one trial for the given step
#

function execute_step {
  echo ${1:0:$2} > "./inputs/${18}.fasta1"

  mkdir models_${2}

  PREVIOUS=0
  MODEL=0
  TRY=0
  let PREVIOUS=${2}-1
  let TRY=MODELS_PER_STEP*${3}
  until [ $MODEL -eq $MODELS_PER_STEP ]
  do
    if [ ${15} -eq 0 ]
    then
      cp models_${PREVIOUS}/model_${MODEL}.pdb model_c.pdb
    fi
    {
      echo_options ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12} ${13} ${14} ${15} ${16} $2 ${18} ${17}
    } > flags

    ../../stud_src/bin/compbio_app.linuxgccrelease @flags -database ../../database
    cp models/model_0.pdb models_${2}/model_${TRY}.pdb
    cp models/score.fsc models_${2}/model_${TRY}.fsc
    let TRY=TRY+1
    let MODEL=MODEL+1
  done
}

#
# Arguments
#
# 1 the current amino acid size
# 2 number of models desired
#
# Selects the $2 best models
#

function select_models {
  STEP=$1
  MODELS_PER_STEP=$2
  
  : > sort.txt
  
  for scoref in models_${STEP}/*.fsc
  do
    sed -n '3 s/[^0-9-]*\(-*[0-9][0-9]*\.*[0-9]*\).*/\1/p' $scoref >> sort.txt
  done
  
  sort -g sort.txt > res.txt

  rm sort.txt

  LAST_SCORE=`sed -n "${MODELS_PER_STEP} s/[^0-9-]*\(-*[0-9][0-9]*\.*[0-9]*\).*/\1/p" res.txt`

  rm res.txt

  ind=0
  while [ -f models_${STEP}/model_${ind}.pdb ]
  do
    CURR_SCORE=`sed -n '3 s/[^0-9-]*\(-*[0-9][0-9]*\.*[0-9]*\).*/\1/p' models_${STEP}/model_${ind}.fsc`
    if [ `echo "$CURR_SCORE > $LAST_SCORE" | bc` -eq 1 ]
    then
      rm models_${STEP}/model_${ind}.pdb
      rm models_${STEP}/model_${ind}.fsc
    fi
    let ind=ind+1
  done

  ind=0
  for modelf in models_${STEP}/model*.pdb
  do
    mv $modelf models_${STEP}/temp${ind}.pdb
    let ind=ind+1
  done

  ind=0
  for modelf in models_${STEP}/temp*.pdb
  do
    mv $modelf models_${STEP}/model_${ind}.pdb
    let ind=ind+1
  done
}

#
# Checking the arguments
#

if [ $# -ne 1 ] || [ ! -d `dirname $0`/$1 ]
then
  echo "Usage: `basename $0` target."
  exit $E_BADARGS
fi  

cd `dirname $0`/$1

#
# Compiling in case of any changes
#

cd ../../stud_src/

#scons bin mode=release

#
# Starting tests
#

cd ../targets/$1

#
# Parameters
#

MIN_STRUC_SIZE=20

MODELS_PER_STEP=4

#
# Read the protein structure
#

SEQUENCE=`sed 's/^>.*//'  "./inputs/${1}.fasta" | grep -v '^$'`

SIZE=${#SEQUENCE}

STEP=$MIN_STRUC_SIZE

#
# First run
#

execute_step $SEQUENCE $STEP 0 2 0 0 0 0 1 2000 2000 2000 4000 50000 1 1 1 $1
execute_step $SEQUENCE $STEP 1 1.5 0 0 0 0 1 2000 2000 2000 4000 50000 1 1 1 $1
execute_step $SEQUENCE $STEP 2 2.5 0 0 0 0 1 2000 2000 2000 4000 50000 1 1 1 $1
execute_step $SEQUENCE $STEP 3 3 0 0 0 0 1 2000 2000 2000 4000 50000 1 1 1 $1

select_models $STEP $MODELS_PER_STEP

#
# Intermediate runs
#

let STEP=STEP+1
until [ $STEP -eq $SIZE ]
do
  execute_step $SEQUENCE $STEP 0 0.5 1 0 0 0 1 2000 2000 2000 4000 50000 0 0 1 $1
  execute_step $SEQUENCE $STEP 1 1 1 0 0 0 1 2000 2000 2000 4000 50000 0 0 1 $1

  select_models $STEP $MODELS_PER_STEP

  let STEP=STEP+1
done 

#
# Last run
#

execute_step $SEQUENCE $STEP 0 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 1 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 2 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 3 1 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 4 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 5 1 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 6 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 7 1 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 8 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 9 1 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 10 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 11 1 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 12 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 13 1 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 14 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 15 1 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 16 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 17 1 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 18 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 19 1 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 20 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 21 1 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 22 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 23 1 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1
execute_step $SEQUENCE $STEP 24 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 1 $1

exit 0

