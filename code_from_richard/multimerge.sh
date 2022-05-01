#!/bin/bash
#
# Merges a sequence of cluster files, writes result to stdout
#
prev=$(mktemp /tmp/multimerge.tmp.XXXX)

for f in $*; do
  next=$(mktemp /tmp/multimerge.tmp.XXXX)
  ./merge $prev $f > $next
  rm $prev
  prev=$next
done

cat $prev
rm $prev
