#!/bin/bash

MUD=$1
num_beds=$2
rings_per_bed=$3
bedStartRing=$4
overlap=$5	# in number of rings

#file names
NC=$6
P1=$7
P2=$8
P3=$9
P4=${10}

###########################

matlab -nodesktop -r "make_senimg0_cylinder(${NC},${P1},${P2},${P3},${P4},${bedStartRing},${num_beds},${rings_per_bed},${overlap},${MUD}); quit"

