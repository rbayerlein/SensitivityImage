#!/bin/bash

P1=\'/home/rbayerlein/Documents/Projects/20210513_Modified_Sensitivity_Image/sens_data_temp/CTAC_201.sen_img\'
P2=\'/home/rbayerlein/Documents/Projects/20210513_Modified_Sensitivity_Image/sens_data_temp/CTAC_201_mumap_kVp-140_size-256x256x646_vox-2.7344x2.7344x3.img\'
P3=\'/home/rbayerlein/Documents/Projects/20210513_Modified_Sensitivity_Image/sens_data_temp/crys_eff_679x840\'		#crys eff map WITH gaps
P4=\'/home/rbayerlein/Documents/Projects/20210513_Modified_Sensitivity_Image/sens_data_temp/plane_eff_679x679\'	#plane eff map WITH gaps

#matlab -nodesktop -nojvm -r "cd /home/rbayerlein/code/explorer-master/reconstruction/lmrecon_senimg/; make_senimg0(${P1},${P2},${P3},${P4}); quit"
matlab -nodesktop -r "make_senimg0(${P1},${P2},${P3},${P4})"