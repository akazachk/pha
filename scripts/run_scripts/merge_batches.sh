#!/bin/bash
# The first argument is the results directory
# The second argument takes values: "best", "test", or "full"
# This means we are merging pha-best.csv, pha-test.csv, or pha-full.csv from the batches

MASTER_RESULTS_DIR="${PHA_DIR}/results"
SCRIPT_DIR="${PHA_DIR}/scripts/run_scripts"
EXECUTABLE="merge_results.sh"

if [ -z "$1" ]
  then RESULTS_DIR="${MASTER_RESULTS_DIR}/batches/test"
else
  RESULTS_DIR="$1"
fi

TMPNAME="pha-test.csv" # default
if [ "$2" = "full" ]
  then TMPNAME="pha-full.csv"
elif [ "$2" = "best" ]
  then TMPNAME="pha-best.csv"
elif [ "$2" = "test" ]
  then TMPNAME="pha-test.csv"
fi

OUT_DIR="${RESULTS_DIR}"
OUTNAME="${OUT_DIR}/${TMPNAME}"

i=0
for batchname in `ls -d ${RESULTS_DIR}/*/`; do
  i=$(( $i + 1 ))  # maintain line count
  if [ "$2" = "full" ]
    then ${SCRIPT_DIR}/${EXECUTABLE} ${batchname} 0 full
  fi
  if [ $i = 1 ]
    then 
      head -n 2 ${batchname}${TMPNAME} > ${OUTNAME}
  fi
  echo "Copying ${TMPNAME} from ${batchname}"
  tail -n +3 ${batchname}${TMPNAME} | grep DONE >> ${OUTNAME}
done

echo "Copying and sorting merged batch result from ${OUTNAME} to ${MASTER_RESULTS_DIR}/${TMPNAME}"
(head -n 2 ${OUTNAME} && (tail -n +3 ${OUTNAME} | sort)) > ${MASTER_RESULTS_DIR}/${TMPNAME}
cp ${MASTER_RESULTS_DIR}/${TMPNAME} ${OUTNAME}
