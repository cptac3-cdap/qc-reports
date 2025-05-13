#!/bin/sh

# set -x

PROG=`readlink -f "$0"`
DIR=`dirname "$PROG"`

if [ "$1" = "" ]; then
  echo `basename "$PROG"` "- Metrics file *.qcmetrics.tsv or SummaryReports directory required." 1>&2
  exit 1
fi

if [ -d "$1" ]; then
  FULLMETRICS=`ls "$1"/*.qcmetrics.tsv 2>/dev/null`
else
  FULLMETRICS="$1"
fi
if [ ! -f "$FULLMETRICS" ]; then
  echo "Can't find metrics file." 1>&2
  exit 1
fi
FULLMETRICS=`readlink -f "$FULLMETRICS"`
COPYBASE=`basename "$FULLMETRICS" .qcmetrics.tsv`
COPYDIR=`dirname "$FULLMETRICS"`
PWDDIR=`readlink -f "$PWD"`
if [ "$COPYDIR" != "$PWDDIR" ]; then
    for extn in .qcmetrics.tsv .sample.txt .label.txt .mayu.tsv\* .peptides.tsv .summary.tsv \*tmt\*.tsv; do
      for f in "$COPYDIR/$COPYBASE"$extn; do
        g=`basename "$f"`
        if [ ! -f "$g" -a -f "$f" ]; then
          cp "$f" .
        fi
      done
    done
    for f in "$COPYDIR/$COPYBASE".qcmetrics.html; do
      g="$COPYBASE".qcmetrics-orig.html
      if [ ! -f "$g" -a -f "$f" ]; then
        cp "$f" "$g"
      fi
    done
    METRICS="$COPYBASE.qcmetrics.tsv"
else
  METRICS="$COPYBASE.qcmetrics.tsv"
fi

BASE=`basename "$METRICS" .qcmetrics.tsv`
SAMPTXT="$BASE.sample.txt"
SAMPCSV="$BASE.sample.csv"
LABELS="$BASE.label.txt"
REPORT="$BASE.qcmetrics.html"
CPTAC_REPORTS=$HOME/projects/cptac3-cdap/cptac-reports
MAKESAMPLE=$CPTAC_REPORTS/make_sample.py
RSCRIPT="conda run -n renv Rscript"

function checkfile {
   if [ ! -f "$1" ]; then
      echo "File $1 required." 1>&2
      exit 1
   fi
}

checkfile "$METRICS"
checkfile "$SAMPTXT"
checkfile "$LABELS"
checkfile "$BASE.summary.tsv"
checkfile "$BASE.peptides.tsv"
checkfile "$BASE"*.tmt*.tsv

if [ -f "$BASE.mayu.tsv.1" ]; then
  cp -f "$BASE.mayu.tsv.1" "$BASE.mayu.tsv"
fi

checkfile "$BASE.mayu.tsv"

rm -f debug.log
if [ ! -s "$SAMPCSV" -o "$SAMPTXT" -nt "$SAMPCSV" ]; then
  awk -F '\t' 'NR > 1 {print $1".mzML.gz"}' "$METRICS" | python3 "$MAKESAMPLE" "$SAMPTXT" "$LABELS" > "$SAMPCSV"
fi
if [ ! -s "$REPORT" -o "$METRICS" -nt "$REPORT" -o "$SAMPTXT" -nt "$REPORT" ]; then
  $RSCRIPT $DIR/generate.R "$METRICS" "$REPORT"
fi
# rm -f "$SAMPCSV"
rm -f *.Rmd

