#!/bin/bash
# rbayerlein August 2021
# script to invoke the sensitivity image creation process for all MUDs until MUD_max which is calculated from the parameter rings_per_bed.

num_beds=6
rings_per_bed=78
bedStartRing=192
overlap=34	# in number of rings
MUD_max=4

# calculate MUD_max
if [[ $rings_per_bed -lt 86 ]]; then
	echo "bed positions could span over no more than 2 units -> MUD_max = 1"
	MUD_max=1
elif [[ $rings_per_bed -gt 85 ]] && [[ $rings_per_bed -lt 170 ]]; then
	echo "bed positions could span over 3 units -> MUD_max = 2"
	MUD_max=2
elif [[ $rings_per_bed -gt 169 ]] && [[ $rings_per_bed -lt 254 ]]; then
	echo "bed positions could span over 4 units -> MUD_max = 3"
	MUD_max=3
else
	echo "bed positions could span over 4 units -> MUD_max = 4"
	MUD_max=4
fi

###########################

#norm coeff:
NC=\'/media/rbayerlein/SSD_09_Reimund/20210827/Multi-Bed_Phantom_Multi-Bed_Phantom_154523/PET/RawData/1.2.156.112605.159303471608576.210827224523.9.6628.91186/1.2.156.112605.159303471608576.210827225240.9.12756.14170.1.nc\'

#output RAW file name of the final combined sensitivity image
P1=\'/media/rbayerlein/data/recon_data/20210827/Multi-Bed_Phantom_Multi-Bed_Phantom_154523/sen_img/CTAC_201.sen_img\'		#output RAW file name of the final combined sensitivity image
#P1=\'/home/rbayerlein/Code/Recon/senimg/dir_temp/CTAC_120_MIN_201.sen_img\'

# Path to CT image 
P2=\'/media/rbayerlein/SSD_09_Reimund/20210827/Multi-Bed_Phantom_Multi-Bed_Phantom_154523/Image/CTAC_201\'
	
# alternative: give path to mu map; if exists: set boolean 'mu_map_exists' to 1
#P2=\'/media/rbayerlein/data/recon_data/20210714/NERVO_ATILLIO_7696939_144717_Raw/sen_img/CTAC_120_MIN_201_mumap_kVp-140_size-256x256x828_vox-2.7344x2.7344x2.344.img\'
#P2=\'/media/rbayerlein/data/recon_data/20200124/sen_img/CTAC_201_mumap_kVp-140_size-256x256x646_vox-2.7344x2.7344x3.img\' 
#P2=\'/media/rbayerlein/data/recon_data/Blank_scan/CTAC_201_mumap_ZEROS_size-256x256x646_vox-2.7344x2.7344x3.img\'

mu_map_exists=0
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
	 sleep 1; ./run_sens.sh ${i} ${num_beds} ${rings_per_bed} ${bedStartRing} ${overlap} ${NC} ${P1} ${P2} ${mu_map_exists} &
done
wait

#combine all bed positions
matlab -nodesktop -r "combine_sen_img_mod_plane_eff(${P1}, ${num_beds}, ${rings_per_bed}, ${bedStartRing}, ${overlap}, ${MUD_max}); quit"
