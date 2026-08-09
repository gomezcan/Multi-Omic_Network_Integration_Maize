#!/bin/bash

for i in *.narrowPeak ; do
# Distance from gene to peak
# save peak with gene targets and position respedt peak:
# If gene in + and dist +, then peak in genebody
# if gene in + and dist -, then peak in gene promoter
# if gene in - and dist +, then peak in genebody
# if gene in - and dist -, then peak in gene promoters 
# Column 17 id dist between TSS and Peak, if the (dis)< 5001. Then we assume 
# its as a targets. 
# -b file, TSS bed like file: chr start end geneid	.	strand

#bedtools closest -a $i -b $1 -D "b" | awk -v OFS='\t' '{ print $11,$12,$13,$14,$15,$16,$17,$4,$1"_"$2"_"$3}' > Dis2TSS.${i//.GEM_events.narrowPeak/.txt};

Name=${i//.GEM_events.narrowPeak/};
Name=${Name//Q30./};

#Name=${Name//_deM_R1/};
# Name=${Name//_deM_R2/};
# Name=${Name//_deM/};


cat $i | awk -v OFS='\t' '{s=(($3-$2)/2); print $1, int($2+s), int($2+s)}' | bedtools closest -a - -b TSS_B73_RefGen_v4.46.bed -D "b" | awk -v OFS='\t' -v tf=$Name '{ print tf,$1":"$2,$7,$9,$10}' > Dis2TSS.${Name}.txt;

done;

cat Dis2TSS*.txt > Net.Dis2TSS.txt;
