#!/bin/sh

for entry in "examples/clustering/input/dbCAN3_new"/*
do
  SUB="GH"
  if [[ "$entry" == *"$SUB"* ]];then
      FN=$(echo "$entry" | rev | cut -d "/" -f 1 | rev)
      FN=$(echo "$FN" | cut -d "." -f 1)
      echo "$FN" 
      #python hmmscan_combine.py -f $FN
      #python hmmscan_nocombine.py -f $FN
      #python hmm_maker.py -f $FN
      python EZ_analysis.py -f $FN
  fi
done
