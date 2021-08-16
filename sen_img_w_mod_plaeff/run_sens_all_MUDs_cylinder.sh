#!/bin/bash
# rbayerlein August 2021
# script to invoke the sensitivity image creation process for all MUDs until MUD_max which is calculated from the parameter rings_per_bed.

num_beds=14
rings_per_bed=84
bedStartRing=20
overlap=42	# in number of rings
MUD_max=4

# calculate MUD_max
if [[ $rings_per_bed -gt 85 ]]; then
	echo "bed positions could span over 3 units -> MUD_max = 2"
	MUD_max=2
elif [[ $rings_per_bed -gt 169 ]]; then
	echo "bed positions could span over 4 units -> MUD_max = 3"
	MUD_max=3
elif [[ $rings_per_bed -gt 253 ]]; then
	echo "bed positions could span over 4 units -> MUD_max = 4"
	MUD_max=4
else
	echo "bed positions could span over no more than 2 units -> MUD_max = 1"
	MUD_max=1
fi

###########################

NC=\'/home/rbayerlein/Code/Recon/senimg/dir_temp/cylinder.1.nc\'
P1=\'/home/rbayerlein/Code/Recon/senimg/dir_temp/CTAC_201.sen_img\'		#output name of the final combined sensitivity image
P2=\'/media/rbayerlein/data/recon_data/20200124/sen_img/CTAC_201_mumap_kVp-140_size-256x256x646_vox-2.7344x2.7344x3.img\' #CT image
#P2=\'/media/rbayerlein/data/recon_data/Blank_scan/CTAC_201_mumap_ZEROS_size-256x256x646_vox-2.7344x2.7344x3.img\'
P3=\'/home/rbayerlein/Code/Recon/senimg/dir_temp/crys_eff_679x840\'		#crys eff map WITH gaps, file WILL BE CREATED at that location
P4=\'/home/rbayerlein/Code/Recon/senimg/dir_temp/plane_eff_679x679\'		#plane eff map WITH gaps, file WILL BE CREATED at that location

###########################

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

for ((i = 0 ; i <= ${MUD_max} ; i++))
do
	./run_sens_cylinder.sh ${i} ${num_beds} ${rings_per_bed} ${bedStartRing} ${overlap} ${NC} ${P1} ${P2} ${P3} ${P4} & 
done

#combine all bed positions
matlab -nodesktop -r "combine_sen_img_mod_plane_eff(${P1}, ${num_beds}, ${rings_per_bed}, ${bedStartRing}, ${overlap}, 1); quit"
