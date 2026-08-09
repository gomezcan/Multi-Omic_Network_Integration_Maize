#!/bin/bash 

### Parametera
## print_bound_seqs: print actual bound sequence. Which is usefull to further 
## motif analysis
##
## --ls: sort ouput by location instead of P val
## --ex: File that contains subset of genome regions to be excluded. Each line contains a region in "chr:start-end" format. NOT IMPLEMENTED YED.
## --fold: Fold (IP/Control) cutoff to filter predicted events (default=3) 
## --min: minimum number of events called by GPS/GEM to learn a read distribution and continue the next round (default=500). 
## 	  Note: with less number of binding events, the read distribution might be more noisy
## --k_seqs: the number of top ranking events to get sequences for motif discovery (default=5000)
## --pp_nmotifs: the max number of top ranking motifs to set the motif-based positional prior (default=1) 
## 

#shopt -s expand_aliases
gem='java -jar -Xmx20G /maindisk/fabio/SoftwareInstalled/gemPC/gem.jar'

SIZES=../GenomeSize_B73.sizes
Control=../DeDup.All.gDNA.B73.bed
dis=../Read_Distribution_default.txt

# Mode 1  control gDNA and relax
echo ".... start Model 1 ..... "

for i in DeDup.Q30.*.bam; do 
out=${i//DeDup./};
out=${out//.bam/};

#$gem --d $dis t --genome ../ChrsB73/ --expt $i --ctrl $Control --t 50 --f BED --g $SIZES --k_min 6 --k_max 15 \
#	--relax --k_seqs 1000 --outNP --sl --out M1_$out;

#echo ".... Done M1 $i ..... "
echo ""
done;


# Mode 2 No control and relax
echo ".... start Model 2 ..... "
for i in DeDup.*.bam; do
out=${i//DeDup.Q30.B73./};
out=${out//.bam/};

$gem --d $dis t --genome ../ChrsB73/ --expt $i --t 50 --f BAM --g $SIZES --k_min 6 --k_max 15 \
	--relax --outNP --sl --out M2_$out;

echo ".... Done M2 $i ..... "
echo ""
done;

#################################    IMPORTANT NOTE    ######################################
# Due to the large variability of datasets, the refined empirical spatial distribution may  #
# become too noisy because too few events are used to estimate the distribution. Running    # 
# GEM further may make the empirical distribution even worse. Therefore, GEM save the       # 
# output files from each round, and allow the user to check the spatial distributions to    #
# decide which one is the best. To facilitate this process, GEM outputs an image to plot    #
# all the read distribution curves (xxx_All_Read_Distributions.png). Ideally, the read      # 
# distribution of later rounds should be smooth and similar to that of round 0. If so, user #
# can use the event predictions from these rounds. If that is not the case, we would        #
# suggest to use the round 0 prediction results                                             #
#############################################################################################

