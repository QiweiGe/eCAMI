#!/bin/sh

for entry in "examples/clustering/input/dbCAN3_new"/*
do
  #specify the CAZy family name by this line
  SUB="GH"
  if [[ "$entry" == *"$SUB"* ]];then
      FN=$(echo "$entry" | rev | cut -d "/" -f 1 | rev)
      FN=$(echo "$FN" | cut -d "." -f 1)
      echo "$FN" 
      
      #extract each protein's domain. The output is in hmm_analysis_combinded. The cut protein is in cut_domain_seq folder. 
      python hmmscan_combine.py -f $FN
      
      #extract each protein's domain. The output is in hmm_analysis_no_combinded. The cut protein is in cut_domain_seq folder.
      #python hmmscan_nocombine.py -f $FN
      
      #build HMM based on protein with domain regions only. The output is in hmm_analysis_combinded or hmm_analysis_no_combinded. The hmm file is in hmm folder.
      python hmm_maker.py -f $FN
      
      #combine all HMM files which belong to the same CAZy family together. The output is in hmm_refe_combine 
      python EZ_analysis.py -f $FN
  fi
done
