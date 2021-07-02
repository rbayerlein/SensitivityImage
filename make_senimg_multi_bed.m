function make_senimg_multi_bed(senimg_name, att_image_name, cryseff_name, plaeff_name, bedStartRing, rings_per_bed, num_beds, overlap_percent)

overlap = overlap_percent/100
lastBedEndRing = bedStartRing + (1+ (num_beds-1)*(1-overlap))*rings_per_bed-1
if (bedStartRing < 1 || bedStartRing > 672 || rings_per_bed < bedStartRing || rings_per_bed > 672 || lastBedEndRing > 672)
    disp('Invalid ring numbers. Cannot define bed position.');
end

senimg_name = senimg_name(1:strfind(senimg_name, '.sen_img')); 
senimg_name = [senimg_name, 'From_', num2str(bedStartRing), '_', num2str(rings_per_bed), 'rings_', num2str(num_beds), 'beds_', num2str(overlap_percent), '_%OL.sen_img']; 

senimg_name
att_image_name
cryseff_name
plaeff_name

num_voxels = 239;   
num_slices = 679;  
imagevoxelsize = 2.85; 

image_size = [num_voxels, num_voxels, num_slices];
voxel_size = [imagevoxelsize, imagevoxelsize, imagevoxelsize];


str1 = '_size-'; 
str2 = '_vox-'; 

ind_ss1 = strfind(att_image_name, str1); 
ind_ss2 = strfind(att_image_name, str2); 
ind_end = strfind(att_image_name, '.img'); 

if isempty(ind_ss1) || isempty(ind_ss2)
	disp('Invalid mu map. Return.'); 
	return; 
else
    disp('Valid mu map found. continue.');
end

fid_attn = fopen(att_image_name, 'rb'); 
att_image = fread(fid_attn, inf, 'float'); 
fclose(fid_attn);

size_str = att_image_name((ind_ss1+length(str1)):(ind_ss2-1));
vox_str = att_image_name((ind_ss2+length(str2)):(ind_end-1));

ind_x1 = strfind(size_str, 'x');
ind_x2 = strfind(vox_str, 'x');

att_image_size = [str2double(size_str(1:(ind_x1(1)-1))) str2double(size_str((ind_x1(1)+1):(ind_x1(2)-1))) str2double(size_str((ind_x1(2)+1):end))]

att_voxel_size = [str2double(vox_str(1:(ind_x2(1)-1))) str2double(vox_str((ind_x2(1)+1):(ind_x2(2)-1))) str2double(vox_str((ind_x2(2)+1):end))]

fid_cryseff = fopen(cryseff_name, 'rb'); 
cryseff_wgap = fread(fid_cryseff, inf, 'float'); 
cryseff_wgap = reshape(cryseff_wgap, 679, 840); 
fclose(fid_cryseff); 

fid_plaeff = fopen(plaeff_name, 'rb'); 
plaeff_wgap = fread(fid_plaeff, inf, 'float'); 
plaeff_wgap = reshape(plaeff_wgap, 679, 679); 
fclose(fid_plaeff); 


padd=genpath('./PETsystem');
addpath(padd);

scanner=buildPET('explorer2000mm_unitedimaging');
obj = scanner;

num_bins_axial_reduced = 84;

num_blockrings = 8;
num_gap = 7;
num_axialxtals = 672;
num_transxtals = 840;
num_axialxtals_wgap = num_axialxtals+num_gap;


fprintf('calculate the sensitivity image:\n');

numpix_attn = 256; %255;
numslices_attn = 828; %672;

numpix_attn = size(att_image, 1); 
numslices_attn = size(att_image, 3); 

%bp = scanner.cal_senimg_uih(plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, blockringdiff_no, num_blockrings, num_bins_axial_reduced);  

bp = scanner.cal_senimg_multi_bed(plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, num_blockrings, num_bins_axial_reduced, bedStartRing, rings_per_bed, num_beds, overlap);
fwrite(fopen(senimg_name, 'w'), bp, 'single');
fclose('all');

ss = ['done bed position with rings ', num2str(bedStartRing), ' to ', num2str(lastBedEndRing), ' with ', num2str(num_beds), ' beds and ', num2str(rings_per_bed), ' rings each']; 
disp(ss); 

%% display sens img

dim_x = 239;
dim_y = dim_x;
dim_z = 679;

slice = reshape(bp(:,round(dim_y/2),:), [dim_x, dim_z]);
imshow(slice, []);
colorbar;
pause(5.0); %wait so that user can quickly look at sensitivity map before continuing.









