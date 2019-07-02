#!/bin/bash
# Arguments: results directory, batch mode flag (0 or 1), run type (best, test, or full)
# All are optional
# Default is ${PHA_DIR}/results/test 0 test
#
# The first argument can be the (full) directory containing the results or batch folders containing results
# Alternatively, it can be a 1 (representing batch mode)
# If the first argument is not 1, then the second argument may be a 1 indicating batch mode
# The argument after the 1 only has an effect if it is "best" or "test", meaning the results are contained in ${CUT_TYPE}-best.csv or ${CUT_TYPE}-test.csv

MASTER_RESULTS_DIR="${PHA_DIR}/results"
RUN_TYPE_STUB="-test"
SCRIPT_DIR="${PHA_DIR}/scripts/run_scripts"
export CUT_TYPE="pha"

# Set the RESULTS_DIR and run in batch mode if so specified
if [ -z "$1" ]
then 
  RESULTS_DIR="${MASTER_RESULTS_DIR}/test"
elif [ "$1" = 1 -o "$2" = 1 ]
then 
  echo "Combining results from batches"
  if [ "$1" != "1" ]
  then 
    RESULTS_DIR="$1"
    ${SCRIPT_DIR}/merge_batches.sh "${RESULTS_DIR}" "$3"
  else
    ${SCRIPT_DIR}/merge_batches.sh "${MASTER_RESULTS_DIR}/batches/test" "$2"
  fi
  exit 1
else
  if [ "$1" = 0 ]
  then 
    RESULTS_DIR="${MASTER_RESULTS_DIR}/test"
  elif [ "$1" != "best" -a "$1" != "test" -a "$1" != "full" ]
  then 
    RESULTS_DIR="$1"
  fi
fi

if [ "$1" = "full" -o "$2" = "full" -o "$3" = "full" ]
then 
  RUN_TYPE_STUB="-full"
elif [ "$1" = "best" -o "$2" = "best" -o "$3" = "best" ]
then 
  RUN_TYPE_STUB="-best"
fi

TMPNAME="${CUT_TYPE}${RUN_TYPE_STUB}"
TMPNAME_EXT=".csv"
OUT_DIR="${RESULTS_DIR}"
OUTNAME="${OUT_DIR}/${TMPNAME}${TMPNAME_EXT}"

if [ "${RUN_TYPE_STUB}" = "-full" ]
then 
  echo "Merging results in ${RESULTS_DIR}"
  # Take the first two rows from one of the files, to get the headings
  head -n 2 ${RESULTS_DIR}/${TMPNAME}0${TMPNAME_EXT} > ${OUTNAME}

  # Now for each of the three files, take every line that contains the word DONE
  tail -n +3 ${RESULTS_DIR}/${TMPNAME}0${TMPNAME_EXT} | grep DONE >> ${OUTNAME}
  tail -n +3 ${RESULTS_DIR}/${TMPNAME}1${TMPNAME_EXT} | grep DONE >> ${OUTNAME}
  tail -n +3 ${RESULTS_DIR}/${TMPNAME}2${TMPNAME_EXT} | grep DONE >> ${OUTNAME}
else
  echo "Copying and sorting results from ${OUTNAME} to ${MASTER_RESULTS_DIR}/${TMPNAME}${TMPNAME_EXT}"
  (head -n 2 ${OUTNAME} && (tail -n +3 ${OUTNAME} | sort)) > ${MASTER_RESULTS_DIR}/${TMPNAME}${TMPNAME_EXT}
fi
