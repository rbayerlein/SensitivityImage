function make_senimg0(nc_file, senimg_name, att_image_path, first_bed_ring, num_beds, rings_per_bed, overlap, MUD, mu_map_exists)

%=======================================================================
% ONLY FOR SIMULATION !!
is_simulation=true;
%=======================================================================

    
    senimg_name
    att_image_path
    mu_map_exists = logical(mu_map_exists)
    
    %% get plane and cristal efficiency and save them to sen_img output folder

    if is_simulation % overwrite crys_eff map with ones
        crys_eff_840x672 = ones(672,840);
        plane_eff_672x672 = ones(672, 672);
    else
        disp(nc_file);
        plane_eff_672x672 = get_plaeff(nc_file, num_beds, first_bed_ring, rings_per_bed, overlap);
        crys_eff_840x672 = get_cryseff(nc_file);
    end
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

    %% get mu map

    if mu_map_exists
        disp('mu map exists.');
        att_image_name=att_image_path;
        
        %extract image info from file name
        find_size=strfind(att_image_name, 'size-')+5;
        dim_x = str2num(att_image_name(find_size:find_size+2));
        dim_y = str2num(att_image_name(find_size+4:find_size+6));
        dim_z = str2num(att_image_name(find_size+8:find_size+10));
        att_image_size = [dim_x dim_y dim_z]
        
        find_vox=strfind(att_image_name, 'vox-')+4;
        find_x = strfind(att_image_name(find_vox:length(att_image_name)), 'x');
        
        vx = str2double(att_image_name(find_vox:find_vox+find_x(1)-2));
        vy = str2double(att_image_name(find_vox+find_x(1):find_vox+find_x(2)-2));
        
        str = att_image_name(find_vox+find_x(2):length(att_image_name));
        for it = 1 : length(str)
            if ~isValidSign(str(it))
                found_NaN = it;
                break
            end
        end
        vz = str2double(str(1:found_NaN-1))
        
        att_voxel_size = [vx,vy,vz]
        
        find_kvp=strfind(att_image_name, 'kVp-')+4;
        CT_kVp = str2num(att_image_name(find_kvp:find_kvp+2))
    else
        disp('create mu map');
        [mu_map, pars] = make_mumap(att_image_path);
        att_image_size = pars.img_size;
        att_voxel_size = pars.vox_size;
        CT_kVp = pars.CT_kVp;

        senimg_name_raw = senimg_name_raw(1:length(senimg_name_raw)-1); % cut the '.' at the end
        att_image_name= [senimg_name_raw, '_mumap_kVp-',num2str(CT_kVp), '_size-', num2str(att_image_size(1)),'x',num2str(att_image_size(2)), 'x', num2str(att_image_size(3)), '_vox-',num2str(att_voxel_size(1)), 'x',num2str(att_voxel_size(2)), 'x', num2str(att_voxel_size(3)), '.img']

        fid_musave = fopen(att_image_name, 'wb'); 
        fwrite(fid_musave, mu_map, 'float'); 
        fclose(fid_musave); 


    end
    %% check validity of att map and read it
    str1 = '_size-'; 
    str2 = '_vox-'; 
    ind_ss1 = strfind(att_image_name, str1); 
    ind_ss2 = strfind(att_image_name, str2); 
    ind_end = strfind(att_image_name, '.img'); 

    if isempty(ind_ss1) || isempty(ind_ss2) || isempty(ind_end)
        disp('Invalid mu map. Return.'); 
        return; 
    else
        disp('Valid mu map found. continue.');
    end

    fid_attn = fopen(att_image_name, 'rb'); 
    att_image = fread(fid_attn, inf, 'float'); 
    fclose(fid_attn);

    %% create scanner object and define essential properties
    padd=genpath('../PETsystem');
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


function bool = isValidSign(inputString)
    bool = false;
    switch inputString
        case '0'
            bool = true;
        case '1'
            bool = true;
        case '2'
            bool = true;
        case '3'   
            bool = true;
        case '4'
            bool = true;
        case '5'
            bool = true;
        case '6'
            bool = true;
        case '7'   
            bool = true;
        case '8'
            bool = true;
        case '9'
            bool = true;
        case '.'
            bool = true;
    end
end

