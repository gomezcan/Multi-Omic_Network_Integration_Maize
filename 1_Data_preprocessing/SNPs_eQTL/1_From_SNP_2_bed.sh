#!/bin/bash

for i in *.GEM_events.narrowPeak; do
	#Q30.B73.MYB124_deM.GEM_events.narrowPeak
	name=${i//Q30.B73./};
	name=${i//Q30./};
	name=${name//.GEM_events.narrowPeak/};
	cat $i | awk -v OFS='\t' -v n=$name '{print $1,$2,$3,n}' | sed 's/chr//g' > $name.bed;
done;
