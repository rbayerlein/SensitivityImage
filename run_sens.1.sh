#!/bin/bash

# usage: define start and end ring number of each bed position starting from 1 as the smallest start position and 672 as the largest ring number
# for combining all to one sensitivity image, define start ring, rings per bed and number of overlapping beds

P1=\'/media/rbayerlein/data/recon_data/20200124/sen_img/84rings_AFOV_crys_eff_ones/CTAC_201.sen_img\'		#output name of the final combined sensitivity image
P2=\'/media/rbayerlein/data/recon_data/20200124/sen_img/CTAC_201_mumap_kVp-140_size-256x256x646_vox-2.7344x2.7344x3.img\' #CT image
#P2=\'/media/rbayerlein/data/recon_data/Blank_scan/CTAC_201_mumap_ZEROS_size-256x256x646_vox-2.7344x2.7344x3.img\'
P3=\'/media/rbayerlein/data/recon_data/20200124/sen_img/crys_eff_679x840\'				#crys eff map WITH gaps. file ACTUALLY needs to be at that location
P4=\'/media/rbayerlein/data/recon_data/20200124/sen_img/plane_eff_679x679\'		#plane eff map WITH gaps. file ACTUALLY needs to be at that location

num_beds=15
rings_per_bed=84
bedStartRing=1
overlap=42		# in numbers of rings

#catch invalid overlap
if [ ${overlap} -ge ${rings_per_bed} ]; then
	echo "invalid overlap. abort."
	exit 1
fi

#catch invalid final bed ring (must be smaller or equal 672, i.e. excluding gap)
let b=${rings_per_bed}-${overlap} n=$num_beds-1 c=n*b start=c+${bedStartRing}
let finalRing=${start}+${rings_per_bed}-1
if [ ${finalRing} -gt 672 ]; then
	echo ${finalRing} ": invalid final ring number. abort"
	exit 1
fi

########################### main part

#for all beds: create sen img 
for ((i = 0 ; i < ${num_beds} ; i++))
do
	#let start=${bedStartRing}+${i}*${rings_per_bed}*${overlap}
	let b=${rings_per_bed}-${overlap} c=${i}*b start=c+${bedStartRing}
	let end=${start}+${rings_per_bed}-1
	echo "Current start and end ring numbers: " ${start} ${end}
	matlab -nodesktop -r "make_senimg0(${P1},${P2},${P3},${P4},${start},${end}); quit"
done

#combine all bed positions
matlab -nodesktop -r "combine_sen_img(${P1}, ${num_beds}, ${rings_per_bed}, ${bedStartRing}, ${overlap});"