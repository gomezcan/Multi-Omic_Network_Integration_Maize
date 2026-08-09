#!/bin/bash


echo Sample TFs Condtion Peaks Targets

for i in *.narrowPeak; do

#Q30.B73.WRKY38_Met.GEM_events.narrowPeak     targets.2kb_TSS.Q30.B73.LBD20_deM.txt

Name=${i//.GEM_events.narrowPeak/}
Name=${Name//Q30./};
#
Peak=$(wc -l $i | cut -d' ' -f1);
Targets=$(cat Dis2TSS.${Name}.txt | awk -v OFS='\t' '{if( sqrt($5*$5) <= 2000) print $0 }' | cut -f3 | sort -u | wc -l);
#
Name2=$(echo $Name | sed $'s/_/ /g')

#echo $Name2
#echo $Name $Name2 $Peak $Targets
echo $Name $Name2 $Peak

done;

