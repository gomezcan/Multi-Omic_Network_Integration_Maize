#!/bin/bash

#source activate FabioPython3.8 

#--walk-length
pecanpy --input uniqFullNets_weighted.txt --weighted --workers 50 --dimensions 50 --walk-length 80 --num-walks 10 \
	--directed --output Pecanpy_uFNetsW.Dim40_WL80_nW10.txt --verbose ;
