#!/bin/bash

source activate FabioPython3.8 


pecanpy --input uniqFullNets.10_11_2021.txt --workers 50 --walk-length 180 --num-walks 10  --directed --output pecanpy_ALL.txt --verbose ;
