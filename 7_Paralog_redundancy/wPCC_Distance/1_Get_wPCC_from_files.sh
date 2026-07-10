#!/bin/bash


Call_wpccDist() {
  # read gene 1 and 2 from pair
  g1=$(echo $1 | cut -d'_' -f1);
  g2=$(echo $1 | cut -d'_' -f2);
 
  # list of file nets for each gene
  grep $g1 filesInwPCCnet > tem.$g1;
  grep $g2 filesInwPCCnet > tem.$g2;

  cat tem.$g1 tem.$g2 > tem12.${g1}.${g2};
  rm tem.$g1;
  rm tem.$g2;
	
  # grep files for genes of interes
  while read line; do
	name=$(echo $line | cut -d'/' -f2 |  sed 's/wPCC.//g' | sed 's/.txt//g');
	val1=$(grep $g1 $line | grep $g2 -);
	val2=$(grep $g2 $line | grep $g1 -);

	echo $name$' '$val1 ;
	echo $name$' '$val2 ;

	done < tem12.${g1}.${g2} > wPCC_${g1}_${g2}.txt;

  rm tem12.${g1}.${g2};
		
  echo " ..  Donde $g1 $g2 .."
}

export -f Call_wpccDist

# Call downloading function
parallel -j 40 Call_wpccDist :::: $1;

#./Call_wpccDist wPCC_net_only_TFs/*.txt; do echo $i; done | grep 'Zm00001d053178' > tem 
