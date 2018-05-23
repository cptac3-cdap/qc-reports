#!/bin/sh

DIR=`dirname $0`

for a in "$@"; do
  b=`basename $a .tsv`
  c="$b.html"
  d="$b.pdf"
  Rscript $DIR/generate.R "$a" "$c"
  Rscript $DIR/generate.R "$a" "$d"
done
