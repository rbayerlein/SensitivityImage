#!/bin/bash

# usage: pass MUD as argument for this sensitivity map calculation
MUD=$1

NC=\'/home/rbayerlein/Code/Recon/senimg/dir_temp/cylinder.1.nc\'
P1=\'/home/rbayerlein/Code/Recon/senimg/dir_temp/CTAC_201.sen_img\'		#output name of the final combined sensitivity image
P2=\'/media/rbayerlein/data/recon_data/20200124/sen_img/CTAC_201_mumap_kVp-140_size-256x256x646_vox-2.7344x2.7344x3.img\' #CT image
#P2=\'/media/rbayerlein/data/recon_data/Blank_scan/CTAC_201_mumap_ZEROS_size-256x256x646_vox-2.7344x2.7344x3.img\'
P3=\'/home/rbayerlein/Code/Recon/senimg/dir_temp/crys_eff_679x840\'		#crys eff map WITH gaps, file WILL BE CREATED at that location
P4=\'/home/rbayerlein/Code/Recon/senimg/dir_temp/plane_eff_679x679\'		#plane eff map WITH gaps, file WILL BE CREATED at that location


num_beds=7
rings_per_bed=84
bedStartRing=1
overlap=50	# in percent

###########################

matlab -nodesktop -r "make_senimg0(${NC},${P1},${P2},${P3},${P4},${bedStartRing},${num_beds},${rings_per_bed},${overlap},${MUD}); quit"

