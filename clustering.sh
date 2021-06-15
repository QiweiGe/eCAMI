#!/bin/sh

for entry in "examples/clustering/input/dbCAN3_new"/*
do
  echo "$entry"
  ##specify CAZy family name by this line
  SUB="GH71"
  if [[ "$entry" == *"$SUB"* ]];then
      FN=$(echo "$entry" | rev | cut -d "/" -f 1 | rev)
      FN=$(echo "$FN" | cut -d "." -f 1)
      python clustering.py -jobs 16 -input $entry -output_dir "examples/clustering/output/dbCAN3_new/"$FN
  fi
done

