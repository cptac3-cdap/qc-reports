#!/bin/sh

DIR=`dirname $0`

if [ "$1" = "--noexpdesign" ]; then
  NOEXPDES="noexpdesign"
  shift;
fi

if [ "$1" = "--dochecks" ]; then
  DOCHECKS="dochecks"
  shift;
fi

for a in "$@"; do
  echo "Processing $a..."
  a=`readlink -f "$a"`
  e=`dirname $a`
  b=`basename $a .tsv`
  c="$e/$b.html"
  d="$e/$b.pdf"
  Rscript $DIR/generate.R "$a" "$c" $NOEXPDES $DOCHECKS || exit 1
  # Rscript $DIR/generate.R "$a" "$d"
done
