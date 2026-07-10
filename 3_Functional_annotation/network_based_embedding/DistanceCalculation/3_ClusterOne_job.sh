#!/bin/bash

#
ClusterONE='java -jar /maindisk/fabio/SoftwareInstalled/cluster_one-1.0.jar'

$ClusterONE -d 0.1 --fluff -s 30 $1 > $2;
