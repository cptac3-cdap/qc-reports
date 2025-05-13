#!/bin/sh

DIR=`dirname "$0"`

function runreport {
  echo "RUNNING: $1"
  $DIR/runadhocreport.sh "$1" && echo "DONE: $1" || echo "ERROR: $1"
}

function allqcmetrics {
  for b in PDC CPTAC; do
    find /data3/$b-CDAP -name "*.qcmetrics.tsv" | fgrep -v attic | fgrep -v forportal | fgrep -v testing | fgrep -v Rename | fgrep -v "SummaryReports.old" | fgrep -v Tutorial | fgrep -v 'CPTAC3_non-ccRCC_JHU_Phosphoproteome-1.qcmetrics.tsv' | fgrep -v 'CPTAC3_UCEC_Confirmatory_Proteome.qcmetrics.tsv' | fgrep 'SummaryReports'
  done
}

for f in `allqcmetrics | sed 's/\/\([^/]*.qcmetrics.tsv\)/:\1/' | sort -t: -k2,2 | tr : /`; do
  runreport "$f"
done

# for f in `find . -name "*.qcmetrics.tsv"`; do
#   runreport "$d"
# done

for f in ./*.qcmetrics.tsv; do
  BASE=`basename "$f" .qcmetrics.tsv`
  if [ ! -f "$BASE.qcmetrics.html" ]; then
    rm -f "$BASE."*
  fi
done

$HOME/bin/generate_directory_index_caddystyle.py --filter "*.html" "*.expdesign.tsv"
