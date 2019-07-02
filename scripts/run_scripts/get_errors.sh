#!/bin/bash
# Arguments: results directory, batch mode flag (0 or 1), run type (best, test, or full)
# All are optional
# Default is ${PHA_DIR}/results/instances/500x500 0 full
#
# The first argument can be the (full) directory containing the results or batch folders containing results
# Alternatively, it can be a 1 (representing batch mode)
# If the first argument is not 1, then the second argument may be a 1 indicating batch mode
# The argument after the 1 only has an effect if it is "best" or "test", meaning the results are contained in ${CUT_TYPE}-best.csv or ${CUT_TYPE}-test.csv

MASTER_RESULTS_DIR="${PHA_DIR}/results/instances"
DEFAULT_RESULTS_DIR="${MASTER_RESULTS_DIR}/500x500"
RUN_TYPE_STUB="-full"
SCRIPT_DIR="${PHA_DIR}/scripts/run_scripts"

if [ -z "${CUT_TYPE}" ]
then 
  echo "*** ERROR: Need to specify cut type."
  exit 1
fi

# Set the RESULTS_DIR and run in batch mode if so specified
if [ -z "$1" ]
  then RESULTS_DIR="${DEFAULT_RESULTS_DIR}"
elif [ "$1" = 1 -o "$2" = 1 ]
  then 
    echo "Combining errors from batches"
    if [ "$1" != "1" ]
      then 
        RESULTS_DIR="$1"
        ${SCRIPT_DIR}/merge_batches_errors.sh "${RESULTS_DIR}" "$3"
    else
      ${SCRIPT_DIR}/merge_batches_errors.sh "${MASTER_RESULTS_DIR}/batches/500x500" "$2"
    fi
    exit 1
else
  if [ "$1" = 0 ]
    then RESULTS_DIR="${DEFAULT_RESULTS_DIR}"
  elif [ "$1" != "best" -a "$1" != "test" ]
    then RESULTS_DIR="$1"
  fi
fi

if [ "$1" = "test" -o "$2" = "test" -o "$3" = "test" ]
  then 
    RUN_TYPE_STUB="-test"
elif [ "$1" = "best" -o "$2" = "best" -o "$3" = "best" ]
  then 
    RUN_TYPE_STUB="-best"

TMPNAME="${RESULTS_DIR}/${CUT_TYPE}"
TMPNAME_EXT=".csv"
OUT_DIR="${RESULTS_DIR}"
OUTNAME="${OUT_DIR}/errors.csv"

echo "Getting errors and sending them to ${OUTNAME}"

if [ "${RUN_TYPE_STUB}" == "-best" -o "${RUN_TYPE_STUB}" == "test" ]
  then
    # Take the first two rows from the file to get the headings
    head -n 2 ${TMPNAME}${RUN_TYPE_STUB}${TMPNAME_EXT} > ${OUTNAME}

    # Now take every line that contains the word ERROR
    tail -n +3 ${TMPNAME}${RUN_TYPE_STUB}${TMPNAME_EXT} | grep ERROR >> ${OUTNAME}
else
  # Take the first two rows from one of the files, to get the headings
  head -n 2 ${TMPNAME}${RUN_TYPE_STUB}0${TMPNAME_EXT} > ${OUTNAME}

  # Now for each of the three files, take every line that contains the word ERROR
  tail -n +3 ${TMPNAME}${RUN_TYPE_STUB}0${TMPNAME_EXT} | grep ERROR >> ${OUTNAME}
  tail -n +3 ${TMPNAME}${RUN_TYPE_STUB}1${TMPNAME_EXT} | grep ERROR >> ${OUTNAME}
  tail -n +3 ${TMPNAME}${RUN_TYPE_STUB}2${TMPNAME_EXT} | grep ERROR >> ${OUTNAME}
fi
