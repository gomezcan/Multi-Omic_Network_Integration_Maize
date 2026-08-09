#!/bin/bash

for i in *.GEM_events.narrowPeak; do 
	name=${i//.GEM_events.narrowPeak/};
	#name=${name//Q30.B73./};
	
	cat $i ${i//_deM./_Met.} | cut -f1-4 | sort -k1,1 -k2,2n | awk -v OFS='\t' '{print $4,$1,$2+50,$3-50,"."}' > $name.saf;
done;
