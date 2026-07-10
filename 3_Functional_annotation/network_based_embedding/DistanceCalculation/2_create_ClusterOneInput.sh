#!/bin/bash

F1="MR_edgesDB_Dim50_WL80_nW10"
F2="MR_edgesDB_Dim50_WL80_nW10_syntenic"
F3="MR_edgesDB_Dim100_WL100_nW50"
F4="MR_edgesDB_Dim100_WL100_nW50_syntenic"

## w50 Dim100 
cat $F3/MR_MI.pecanpy.* | awk -v OFS='\t' '{if($6>=0.01) print $0}' | grep -v 'GeneID' | cut -f1-2,5 > InputClusterONE_Dim100_WL100_nW50.txt;

# w50 Dim100 syntenic
cat $F4/MR_MI.pecanpy.* | awk -v OFS='\t' '{if($6>=0.01) print $0}' | grep -v 'GeneID' | cut -f1-2,5 > InputClusterONE_Dim100_WL100_nW50_syntenic.txt;

## w50 Dim50
cat $F1/MR_MI.pecanpy.* | awk -v OFS='\t' '{if($6>=0.01) print $0}' | grep -v 'GeneID' | cut -f1-2,5 > InputClusterONE_Dim50_WL80_nW10.txt;

## w50 Dim50 syntenic
cat $F2/MR_MI.pecanpy.* | awk -v OFS='\t' '{if($6>=0.01) print $0}' | grep -v 'GeneID' | cut -f1-2,5 > InputClusterONE_Dim50_WL80_nW10_syntenic.txt;

## based on MI
cat $2/MR_MI.pecanpy.* | awk -v OFS='\t' '{if($3>=0.5) print $0}' | grep -v 'GeneID' | cut -f1-3 > MI0.5_Net.txt

## based on MR
cat $F2/MR_MI.pecanpy.* | awk -v OFS='\t' '{if($6>=0.005) print $0}' | grep -v 'GeneID' | cut -f1-2,5 > InputClusterONE_Dim50_WL80_nW10_0.005_syntenic.txt
