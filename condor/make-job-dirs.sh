#! /bin/bash

for ((i=1; i<=$1; i++))
do
  mkdir -p beautydata-$3-$4/$i
  echo "loo\n$2\n$3\n$4\n$i" > beautydata-$3-$4/$i/params.txt
done
