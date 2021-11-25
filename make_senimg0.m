function make_senimg0(senimg_name, att_image_name, cryseff_name, plaeff_name, bedStartRing, bedEndRing)

if (bedStartRing < 1 || bedStartRing > bedEndRing || bedEndRing < bedStartRing || bedEndRing > 672)
    disp('Invalid ring numbers. Cannot define bed position.');
end

senimg_name = senimg_name(1:strfind(senimg_name, '.sen_img')); 
senimg_name = [senimg_name, num2str(bedStartRing), '_', num2str(bedEndRing), '.sen_img']; 

disp(senimg_name);
disp(att_image_name);
disp(cryseff_name);
disp(plaeff_name);

%% PET image parameters
num_voxels = 239;       % x and y 
num_slices = 679;       % z
imagevoxelsize = 2.85;  % isotropical

image_size = [num_voxels, num_voxels, num_slices];
voxel_size = [imagevoxelsize, imagevoxelsize, imagevoxelsize];

%% find mu map image parameters from file name
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

size_str = att_image_name((ind_ss1+length(str1)):(ind_ss2-1));
vox_str = att_image_name((ind_ss2+length(str2)):(ind_end-1));

ind_x1 = strfind(size_str, 'x');
ind_x2 = strfind(vox_str, 'x');

att_image_size = [str2double(size_str(1:(ind_x1(1)-1))) str2double(size_str((ind_x1(1)+1):(ind_x1(2)-1))) str2double(size_str((ind_x1(2)+1):end))]
att_voxel_size = [str2double(vox_str(1:(ind_x2(1)-1))) str2double(vox_str((ind_x2(1)+1):(ind_x2(2)-1))) str2double(vox_str((ind_x2(2)+1):end))]

% check image dimensions and rebin if neccessary:

if (att_image_size(1) ~= image_size(1) || att_voxel_size(1) ~= voxel_size(1) || att_image_size(2) ~= image_size(2) || att_voxel_size(2) ~= voxel_size(2) || att_image_size(3) ~= image_size(3) || att_voxel_size(3) ~= voxel_size(3))
    att_image_name = make_Downsampling(att_image_name);
    att_image_size = image_size;
    att_voxel_size = voxel_size;
end

%% read attn image
fid_attn = fopen(att_image_name, 'rb'); 
att_image = fread(fid_attn, inf, 'float'); 
fclose(fid_attn);


%% read in plane and crystal efficiency

fid_cryseff = fopen(cryseff_name, 'rb'); 
cryseff_wgap = fread(fid_cryseff, inf, 'float'); 
cryseff_wgap = reshape(cryseff_wgap, 679, 840); 
fclose(fid_cryseff); 

fid_plaeff = fopen(plaeff_name, 'rb'); 
plaeff_wgap = fread(fid_plaeff, inf, 'float'); 
plaeff_wgap = reshape(plaeff_wgap, 679, 679); 
fclose(fid_plaeff); 

%% build scanner
padd=genpath('./PETsystem');
addpath(padd);

scanner=buildPET('explorer2000mm_unitedimaging');
obj = scanner;

num_bins_axial_reduced = 84;
num_blockrings = 8;

%% calculate sensitivity image
fprintf('calculate the sensitivity image: \n');
tStart = tic; % start timer
bp = scanner.cal_senimg_single_bed(plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, num_blockrings, num_bins_axial_reduced, bedStartRing, bedEndRing);
fwrite(fopen(senimg_name, 'w'), bp, 'single');
fclose('all');
% display time elapsed since start of calculation
tEnd = toc(tStart);
ss = ['done bed position with rings ', num2str(bedStartRing), ' to ', num2str(bedEndRing)]; 
disp(ss); 
fprintf('===\nTime elapsed for this calculation: %0.0f min, %0.2f sec \n===\n', floor(tEnd/60), rem(tEnd, 60));

%% display sens img

dim_x = num_voxels;
dim_y = dim_x;
dim_z = num_slices;

slice = reshape(bp(:,round(dim_y/2),:), [dim_x, dim_z]);
imshow(slice, []);
colorbar;
pause(5.0); % wait so that user can quickly look at sensitivity map before continuing.
