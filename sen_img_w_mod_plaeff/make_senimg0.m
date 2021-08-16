function make_senimg0(nc_file, senimg_name, att_image_path, first_bed_ring, num_beds, rings_per_bed, overlap, MUD)

    nc_file
    senimg_name
    att_image_path

    %% get plane and cristal efficiency and save them to sen_img output folder
    plane_eff_672x672 = get_plaeff(nc_file, num_beds, first_bed_ring, rings_per_bed, overlap);
    crys_eff_840x672 = get_cryseff(nc_file);
    [crys_eff_840x679, plane_eff_679x679] = add_axial_gap(crys_eff_840x672, plane_eff_672x672); 


    cryseff_name = [senimg_name(1:find_last_slash(senimg_name)), 'crys_eff_679x840']
    fid_cryseff = fopen(cryseff_name, 'wb');
    fwrite(fid_cryseff, crys_eff_840x679, 'float');
    fclose(fid_cryseff);

    plaeff_name = [senimg_name(1:find_last_slash(senimg_name)), 'plane_eff_679x679']
    fid_plaeff = fopen(plaeff_name, 'wb');
    fwrite(fid_plaeff, plane_eff_679x679, 'float');
    fclose(fid_plaeff);

    %% define output name    
    senimg_name_raw = senimg_name(1:strfind(senimg_name, '.sen_img'));
    senimg_name = [senimg_name_raw, num2str(num_beds), 'beds_', num2str(rings_per_bed), 'rings_from', num2str(first_bed_ring), '_', num2str(overlap), 'rings_overlap_MUD', num2str(MUD), '.sen_img']

    last_bed_ring = first_bed_ring + rings_per_bed + (num_beds-1)*(rings_per_bed-overlap)-1
    num_voxels = 239;   
    num_slices = 679;  
    imagevoxelsize = 2.85; 

    image_size = [num_voxels, num_voxels, num_slices];
    voxel_size = [imagevoxelsize, imagevoxelsize, imagevoxelsize];

    %% create mu map

    [mu_map, pars] = make_mumap(att_image_path);
    att_image_size = pars.img_size;
    att_voxel_size = pars.vox_size;
    CT_kVp = pars.CT_kVp;

    senimg_name_raw = senimg_name_raw(1:length(senimg_name_raw)-1); % cut the '.' at the end
    att_image_name= [senimg_name_raw, '_mumap_kVp-',num2str(CT_kVp), '_size-', num2str(att_image_size(1)),'x',num2str(att_image_size(2)), 'x', num2str(att_image_size(3)), '_vox-',num2str(att_voxel_size(1)), 'x',num2str(att_voxel_size(2)), 'x', num2str(att_voxel_size(3)), '.img']

    fid_musave = fopen(att_image_name, 'wb'); 
    fwrite(fid_musave, mu_map, 'float'); 
    fclose(fid_musave); 

    str1 = '_size-'; 
    str2 = '_vox-'; 

    %% check validity of att map and read it
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

    %% create scanner object and define essential properties
    padd=genpath('/home/rbayerlein/Code/Recon/senimg/PETsystem');
    addpath(padd);

    scanner=buildPET('explorer2000mm_unitedimaging');
    obj = scanner;

    num_bins_axial_reduced = 84;
    num_blockrings = 8;

    %% start main process
    fprintf('calculate the sensitivity image:\n');


    bp = scanner.cal_senimg_uih(plane_eff_679x679, crys_eff_840x679, att_image, att_image_size, att_voxel_size, image_size, voxel_size, MUD, num_blockrings, num_bins_axial_reduced, first_bed_ring, last_bed_ring);  

    %bp = scanner.cal_senimg_single_bed(plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, num_blockrings, num_bins_axial_reduced, first_bed_ring, bedEndRing);
    fwrite(fopen(senimg_name, 'w'), bp, 'single');
    fclose('all');

    ss = ['done bed position with MUD ', num2str(MUD)]; 
    disp(ss); 

    %% display sens img

    dim_x = 239;
    dim_y = dim_x;
    dim_z = 679;

    slice = reshape(bp(:,round(dim_y/2),:), [dim_x, dim_z]);
    imshow(slice, []);
    colorbar;
    pause(5.0); %wait so that user can quickly look at sensitivity map before continuing.
    quit;

end

%% function to find last slash in a string
function pos = find_last_slash(inputString)
    pos = 1;
    for k = 1 : length(inputString)
        if strcmp(inputString(k), '/')
            pos = k;
        end
    end
end




