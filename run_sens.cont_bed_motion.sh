#!/bin/bash

# usage: define start and end ring number of each bed position starting from 1 as the smallest start position and 679 as the largest ring number
# for combining all to one sensitivity image, define start ring, rings per bed and number of overlapping beds

P1=\'/mnt/data/rbayerlein/explorer/20200124/cylinder_phantom3_E-20200124-172216-11_172342/UCD/Image/cont_bed_motion/CTAC_201.sen_img\'		#output name of the final combined sensitivity image
P2=\'/mnt/data/rbayerlein/explorer/20200124/cylinder_phantom3_E-20200124-172216-11_172342/UCD/Image/CTAC_201_mumap_kVp-140_size-256x256x646_vox-2.7344x2.7344x3.img\' #CT image
P3=\'/mnt/data/rbayerlein/explorer/20200124/cylinder_phantom3_E-20200124-172216-11_172342/UCD/Image/crys_eff_679x840\'		#crys eff map WITH gaps
P4=\'/mnt/data/rbayerlein/explorer/20200124/cylinder_phantom3_E-20200124-172216-11_172342/UCD/Image/plane_eff_679x679\'		#plane eff map WITH gaps

num_beds=168
rings_per_bed=84
bedStartRing=1

###########################

#for i in 1..${num_beds}
for ((i = 0 ; i < ${num_beds} ; i++))
do
	#let start=${bedStartRing}+${i}*${rings_per_bed}*${overlap}
	let start=${bedStartRing}+${i}
	let end=${start}+${rings_per_bed}-1
	echo "Current start and end ring numbers: " ${start} ${end}
	matlab -nodesktop -r "make_senimg0(${P1},${P2},${P3},${P4},${start},${end}); quit"
done

#start=1
#end=80
#matlab -nodesktop -r "make_senimg0(${P1},${P2},${P3},${P4},${start},${end}); quit"

matlab -nodesktop -r "combine_sen_img_cont_bed_motion(${P1}, ${num_beds}, ${rings_per_bed}, ${bedStartRing});"
