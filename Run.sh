#!/bin/bash

E_BADARGS=65

MIN_STRUC_SIZE=20

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
    echo "-in:file:native model_0.pdb"
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

scons bin mode=release

#
# Starting tests
#

cd ../targets/$1

#
# Read the protein structure
#

SEQUENCE=`sed 's/^>.*//'  "./inputs/${1}.fasta" | grep -v '^$'`

SIZE=${#SEQUENCE}

#
# First run
#

echo ${SEQUENCE:0:$MIN_STRUC_SIZE} > "./inputs/${1}.fasta1"

{
  echo_options 2 0 0 0 0 0 2000 2000 2000 4000 50000 1 1 $MIN_STRUC_SIZE $1 1
} > flags

../../stud_src/bin/compbio_app.linuxgccrelease @flags -database ../../database
cp models/model_0.pdb model_0.pdb

#
# Intermediate runs
#

COUNT=$MIN_STRUC_SIZE
let COUNT=COUNT+1
until [ $COUNT -eq $SIZE ]
do
  echo ${SEQUENCE:0:$COUNT} > "./inputs/${1}.fasta1"

  {
    echo_options 0.5 1 0 0 0 1 200 200 200 400 5000 0 0 $COUNT $1 1
  } > flags

  ../../stud_src/bin/compbio_app.linuxgccrelease @flags -database ../../database
  cp models/model_0.pdb model_0.pdb

  let COUNT=COUNT+1
done 

#
# Last run
#

echo ${SEQUENCE:0:$SIZE} > "./inputs/${1}.fasta1"

{
  echo_options 0.5 1 0 0 0 0 2000 2000 2000 4000 50000 0 0 $SIZE $1 100
} > flags

../../stud_src/bin/compbio_app.linuxgccrelease @flags -database ../../database
cp models/model_0.pdb model_0.pdb


exit 0

