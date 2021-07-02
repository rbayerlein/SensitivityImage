#!/bin/bash

# usage: define start ring, number of rings per bed, number of beds and overlap in percent
# for combining all to one sensitivity image, define start ring, rings per bed and number of overlapping beds

P1=\'/media/rbayerlein/data/recon_data/Blank_scan/CTAC_201.BLANK.sen_img\'		#output name of the final combined sensitivity image
P2=\'/media/rbayerlein/data/recon_data/Blank_scan/CTAC_201_mumap_ZEROS_size-256x256x646_vox-2.7344x2.7344x3.img\' #CT image
P3=\'/media/rbayerlein/data/recon_data/20200124/sen_img/crys_eff_679x840\'		#crys eff map WITH gaps
P4=\'/media/rbayerlein/data/recon_data/20200124/sen_img/plane_eff_679x679\'		#plane eff map WITH gaps

num_beds=2
rings_per_bed=60
bedStartRing=1
overlap=50	# in percent

###########################

matlab -nodesktop -r "make_senimg_multi_bed(${P1},${P2},${P3},${P4},${bedStartRing},${rings_per_bed},${num_beds},${overlap}); quit"

#start=1
#end=80
#matlab -nodesktop -r "make_senimg0(${P1},${P2},${P3},${P4},${start},${end}); quit"

#matlab -nodesktop -r "combine_multi_bed_sen_img(${P1}, ${num_beds}, ${rings_per_bed}, ${bedStartRing}, ${overlap});"