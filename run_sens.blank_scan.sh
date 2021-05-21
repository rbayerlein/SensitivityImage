#!/bin/bash

# usage: define start and end ring number of each bed position starting from 1 as the smallest start position and 679 as the largest ring number
# for combining all to one sensitivity image, define start ring, rings per bed and number of overlapping beds

P1=\'/home/rbayerlein/Documents/Projects/20210513_Modified_Sensitivity_Image/sens_data_temp/CTAC_201.BLANK.sen_img\'		#output name of the final combined sensitivity image
P2=\'/home/rbayerlein/Documents/Projects/20210513_Modified_Sensitivity_Image/sens_data_temp/CTAC_201_mumap_ZEROS_size-256x256x646_vox-2.7344x2.7344x3.img\' #CT image
P3=\'/home/rbayerlein/Documents/Projects/20210513_Modified_Sensitivity_Image/sens_data_temp/crys_eff_679x840\'		#crys eff map WITH gaps
P4=\'/home/rbayerlein/Documents/Projects/20210513_Modified_Sensitivity_Image/sens_data_temp/plane_eff_679x679\'		#plane eff map WITH gaps

start=1
end=80
matlab -nodesktop -r "make_senimg0(${P1},${P2},${P3},${P4},${start},${end}); quit"

num_beds=1
ring_per_bed=80
bedStartRing=1
matlab -nodesktop -r "combine_sen_img(${P1}, ${num_beds}, ${ring_per_bed}, ${bedStartRing});"