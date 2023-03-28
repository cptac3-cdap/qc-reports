#!/bin/sh

if [ $# -lt 3 ]; then
  echo "Usage: qcmetrics.sh <mzMLfile> <mzIdentMLfile> <qcmetricsfile>" 1>&2
  exit 1;
fi

BASE=`readlink -f "$0"`
BASE=`dirname "$BASE"`

METRICR="$BASE/qcmetrics.r"

# For Edwards lab installation of R and java
export PATH="/tools/bin:$PATH"

RSCRIPT=${RSCRIPT:-Rscript} 

ulimit -s unlimited
( "$RSCRIPT" "$METRICR" "$1" "$2" "$3" >/dev/null 2>&1 && test -s "$3" ) || exit 1
