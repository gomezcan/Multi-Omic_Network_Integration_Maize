#!/bin/bash

Full=CoExp_NetworkFinal.10_11_2021.txt

Nets=$(cut -f1 $Full | sort -u | grep -v 'Net')

for i in $Nets; do
	#
	grep -w $i $Full | cut -f2,3 > sub
	echo "ready to $i"
	#
	pecanpy --input sub --workers 40 --directed --output temout --verbose ;
	#
	cat temout | awk '(NR>1)' > CoExp.${i}.pecanpy.out.txt
	rm temout;
	rm sub;

done;

#nohup pecanpy --input PDI.pecanpy.input.txt --workers 30 --directed --output PDI.pecanpy.out.txt & 
