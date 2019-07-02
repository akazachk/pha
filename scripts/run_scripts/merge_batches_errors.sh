#!/bin/bash
# The first argument is the results directory
# The only values for the second argument having any effect are "best" or "test", which will mean we are merging ${CUT_TYPE}-best.csv (${CUT_TYPE}-test.csv) from the batches
# For batches, keep the list of batches in batch_list.txt under scripts/run_scripts, or change the variables below

MASTER_RESULTS_DIR="${PHA_DIR}/results/instances"
BATCHLIST_DIR="${PHA_DIR}/scripts/run_scripts"
BATCHLIST="${BATCHLIST_DIR}/batch_list.txt"
SCRIPT_DIR="${PHA_DIR}/scripts/run_scripts"
EXECUTABLE="get_errors.sh"

if [ -z "${CUT_TYPE}" ]
then 
  echo "***ERROR: Need to specify cut type."
  exit 1
fi

if [ -z "$1" ]
  then RESULTS_DIR="${MASTER_RESULTS_DIR}/batches"
else
  RESULTS_DIR="$1"
fi

TMPNAME="${CUT_TYPE}-full.csv"
if [ "$2" = "best" ]
  then TMPNAME="${CUT_TYPE}-best.csv"
elif [ "$2" = "test" ]
  then TMPNAME="${CUT_TYPE}-test.csv"
fi
OUT_DIR="${RESULTS_DIR}"
OUTNAME="${OUT_DIR}/errors.csv"

i=0
exec <${BATCHLIST}  # redirect file into our stdin
while read ; do    # read each line into REPLY variable
  i=$(( $i + 1 ))  # maintain line count
  batch_name=${REPLY}
  if [ "$2" != "best" -a "$2" != "test" ]
    then ${SCRIPT_DIR}/${EXECUTABLE} ${RESULTS_DIR}/${batch_name}
  fi
  if [ $i = 1 ]
    then head -n 2 ${RESULTS_DIR}/${batch_name}/${TMPNAME} > ${OUTNAME}
  fi
  echo "Copying errors in ${batch_name}/${TMPNAME} to ${OUTNAME}"
  tail -n +3 ${RESULTS_DIR}/${batch_name}/${TMPNAME} | grep ERROR >> ${OUTNAME}
done

