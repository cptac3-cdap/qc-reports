#!/bin/sh

DIR=`dirname $0`

for a in "$@"; do
  a=`readlink -f "$a"`
  e=`dirname $a`
  b=`basename $a .tsv`
  c="$e/$b.html"
  d="$e/$b.pdf"
  Rscript $DIR/generate.R "$a" "$c"
  Rscript $DIR/generate.R "$a" "$d"
done
