% A class to perform fast listmode- or sinogram-based PET image reconstruction
%
% Author: Jian Zhou
% Author: Xuezhu Zhang
%

classdef PETsystem
    
    properties(GetAccess = 'public', SetAccess = 'private')
        name_tag; % scanner name
        system_parms; % structure that contains various system parameters
    end
    
    
    
    methods(Access = 'public')
        
        
        
        
        function obj = PETsystem(...
                tag, ...
                ring_diameter, ...
                crystal_size, ...
                crystal_gap_size, ...
                crystal_array_size, ...
                number_of_detector_modules_transaxial, ...
                number_of_detector_modules_axial, ...
                number_of_DOI_bins, ...
                detector_module_initial_angle_offset, ...
                detector_modula_axial_extra_offsets, ...
                number_of_projections_per_angle, ...
                tof_info)
            
            obj.name_tag = tag;
            obj.system_parms.ring_diameter = ring_diameter;
            obj.system_parms.crystal_size = crystal_size;
            obj.system_parms.crystal_gap_size = crystal_gap_size;
            obj.system_parms.crystal_array_size = crystal_array_size;
            obj.system_parms.number_of_detector_modules_transaxial = ...
                number_of_detector_modules_transaxial;
            obj.system_parms.number_of_detector_modules_axial = ...
                number_of_detector_modules_axial;
            obj.system_parms.number_of_DOI_bins = number_of_DOI_bins;
            obj.system_parms.detector_module_initial_angle_offset = ...
                detector_module_initial_angle_offset;
            obj.system_parms.detector_module_axial_extra_offsets = ...
                detector_modula_axial_extra_offsets;
            obj.system_parms.number_of_projections_per_angle = ...
                number_of_projections_per_angle;
            obj.system_parms.tof_info = tof_info;
            obj.system_parms.projector = 'Siddon';
            obj.system_parms.depth_ratio = 0.5;
        end
        
        
        
        
        
        function obj = setProjector(obj, projector)
            %set projector type: `Siddon'
            %function obj = setProjector(projector)
            obj.system_parms.projector = projector;
        end
        
        
        
        
        
        
        function obj = setDepthRatio(obj, ratio)
            %set a ratio to determine average depth interaction point inside crystal
            %0.5 is default, meaning at the crystal center
            %function obj = setDepthRatio(ratio)
            if ratio < 0 || ratio > 1.0
                error('ratio must be in [0 1]!');
            end
            obj.system_parms.depth_ratio = ratio;
        end
        
        
        
        
        
        function obj = setNumberOfProjectionsPerAngle(obj, num_of_radial_bins)
            %set number of radial bins
            %function obj = setNumberOfProjectionsPerAngle(num_of_radial_bins)
            obj.system_parms.number_of_projections_per_angle = num_of_radial_bins;
        end
        
        
        
        
        
        function dmo = getDetectorModuleAxialOffsets(obj)
            %calculate detector module axial offsets
            %function dmo = getDetectorModuleAxialOffsets
            crystal_axial_pitch = obj.system_parms.crystal_size(3) + ...
                obj.system_parms.crystal_gap_size(3);
            detector_module_axial_size = crystal_axial_pitch * ...
                obj.system_parms.crystal_array_size(2);
            
            nb = obj.system_parms.number_of_detector_modules_axial;
            dmo = (- nb * 0.5 + 0.5 + (0 : (nb-1))) * ...
                detector_module_axial_size + ...
                obj.system_parms.detector_module_axial_extra_offsets;
        end
        
        
        
        
        
        
        function nr = getNumberOfCrystalRings(obj)
            %get number of crystal rings
            %function nr = getNumberOfCrystalRings
            nr = obj.system_parms.crystal_array_size(2) * ...
                obj.system_parms.number_of_detector_modules_axial;
        end
        
        
        
        
        
        
        function na = getDefaultNumberOfAngles(obj)
            %get number of projection angles in default mode
            %function na = getDefaultNumberOfAngles
            na = obj.system_parms.crystal_array_size(1) * ...
                obj.system_parms.number_of_detector_modules_transaxial / 2;
        end
        
        
        
        
        
        
        function ro = getCrystalRingOffsets(obj)   % define axial gap (2014-0311)
            %get crystal ring axial offsets
            %function ro = getCrystalRingOffsets
            crystal_axial_pitch = obj.system_parms.crystal_size(3) + ...
                obj.system_parms.crystal_gap_size(3);
            dm_offsets = getDetectorModuleAxialOffsets(obj);
            nb = obj.system_parms.crystal_array_size(2);
            xt_centers = (- nb * 0.5 + 0.5 + (0 : (nb-1))) * crystal_axial_pitch;
            ro = [];
            for n = 1 : obj.system_parms.number_of_detector_modules_axial
                ro = [ro , xt_centers + dm_offsets(n)];
            end
        end
        
        
        
        
        
        
        function tc = getCrystalTransaxialLocations(obj)
            %calculate crystal bin transaxial coordinates (bin means DOI bin)
            %function tc = getCrystalTransaxialLocations
            tfs = (obj.system_parms.crystal_size(1) + obj.system_parms.crystal_gap_size(1));
            trs = (obj.system_parms.crystal_size(2) + obj.system_parms.crystal_gap_size(2));
            nxtal_trans = obj.system_parms.crystal_array_size(1);
            df = (-nxtal_trans*0.5 + 0.5 + (0 : nxtal_trans-1)) * tfs;
            dr = (-nxtal_trans*0.5 + 0.5 + (0 : nxtal_trans-1)) * trs;
            
            num_of_doi_bins = obj.system_parms.number_of_DOI_bins;
            xtal_loc = zeros(2, num_of_doi_bins, nxtal_trans);
            xtal_size_depth = obj.system_parms.crystal_size(4);
            for n=1: nxtal_trans
                l = df(n)-dr(n);
                dl = l / num_of_doi_bins;
                cl = l * 0.5 - ((0:num_of_doi_bins-1) + 0.5) * dl;
                h = xtal_size_depth;
                dh = h / num_of_doi_bins;
                ch = h * 0.5 - ((0:num_of_doi_bins-1) + 0.5) * dh;
                xtal_loc(:, :, n) = [cl(:) + dr(n), ch(:)]';
            end
            
            R = obj.system_parms.ring_diameter * 0.5;
            nblock_trans = obj.system_parms.number_of_detector_modules_transaxial;
            tc = zeros(2, num_of_doi_bins, nxtal_trans, nblock_trans);
            for n=1:nblock_trans
                % always start at 3 o'clock
                a0 = (n-1) * 2*pi / nblock_trans  + ...
                    obj.system_parms.detector_module_initial_angle_offset * pi / 180;
                xs = xtal_loc(1,:,:);
                ys = xtal_loc(2,:,:);
                
                t = pi * 0.5 + a0;
                x = xs * cos(t) - ys * sin(t) + (R + xtal_size_depth*obj.system_parms.depth_ratio) * cos(a0);
                y = xs * sin(t) + ys * cos(t) + (R + xtal_size_depth*obj.system_parms.depth_ratio) * sin(a0);
                
                tc(1,:,:,n) = x;
                tc(2,:,:,n) = y;
            end
            tc = reshape(tc, 2, num_of_doi_bins, nxtal_trans*nblock_trans);
            tc = squeeze(tc);
        end
        
        
        
        
    end
    
    
    
    
    
    methods
        
        
        
        
        
        function xtal_pairs = getDefaultSinogramCrystalPairs(obj)
            % create crystal pairs according my own numbering scheme
            %function xtal_pairs = getDefaultCrystalPairs
            %
            
            % always equal to total number of crystals per ring divided by
            % 2
            nxtal_trans = obj.system_parms.crystal_array_size(1);   % 35
            xtal_num_trans_total = obj.system_parms.number_of_detector_modules_transaxial * ...
                obj.system_parms.crystal_array_size(1); % 24 * 35 = 840
            
            num_of_angles = xtal_num_trans_total / 2;   % = 420
            
            if obj.system_parms.number_of_projections_per_angle > (fix(xtal_num_trans_total/2)*2 - 2)
                error('too many projections!');
            end
            
            % case when detector module is not located exactly at 3-clock
            if obj.system_parms.detector_module_initial_angle_offset ~= 0
                disp('The 1st detector module is rotated by an angle!');
                disp('NOTE: For crystal pairing, this angle offset is never used directly!');
                disp('but always assume it is equal to 180 / (number_of_detector_modules_transaxial)');
                id0 = 0;
                id1 = num_of_angles;
                nr = fix(num_of_angles/2) * 2 - 2;
                h0 = zeros(nr,1);
                h1 = zeros(nr,1);
                
                for i=0:nr-1
                    
                    id=id0 + fix(nr/2) - i - 1;
                    if id < 0
                        id = id + xtal_num_trans_total;
                    end
                    
                    h0(i+1)=id;
                    id = id1-fix(nr/2) + i;
                    h1(i+1)=id;
                end
            else
                
                id0 = fix(nxtal_trans / 2); % fix(35/2) = 17
                odd = mod(nxtal_trans, 2) ~= 0;
                
                if odd
                    id1 = fix(xtal_num_trans_total / 2) + fix(nxtal_trans / 2); % 420+17 = 437
                    nr = fix(num_of_angles / 2) * 2 - 1; % 419
                else
                    id1 = fix(xtal_num_trans_total / 2) + fix(nxtal_trans / 2) - 1;
                    nr = fix(num_of_angles / 2) * 2;
                end
                
                h0 = zeros(nr,1);   % 419 lines, 1 column
                h1 = zeros(nr,1);
                
                for i=0:nr-1    % loop over all 419 angles and define pairs between opposite crystals h0 and h1
                    if odd
                        id = id0 + fix(nr/2) - i; % 17  + fix(419/2) -i = 226 - i
                    else
                        id = id0 + fix(nr/2) -1 - i;
                    end
                    
                    if id < 0
                        id = id + xtal_num_trans_total;
                    end
                    
                    h0(i+1) = id;
                    if odd
                        id = id1 - fix(nr/2) + i; % 437 - fix(419/2) + i = 228 + i
                    else
                        id = id1 - fix(nr/2) + 1 + i;
                    end
                    h1(i+1) = id;
                end
            end
            
            % pairing
            c = 1; k = 1;
            while c < length(h0) % c < 419
                xtal_id1 = h0(c);
                xtal_id2 = h1(c);
                xp_first_angle(:,k) = [xtal_id1; xtal_id2];
                k = k + 1;
                if (c+1) <= length(h1)
                    xtal_id1 = h0(c);
                    xtal_id2 = h1(c+1);
                    xp_first_angle(:,k) = [xtal_id1; xtal_id2];
                    k = k + 1;
                end
                c = c + 1;
            end
            
            %
            nn = size(xp_first_angle,2);
            num_of_projs_per_angle = obj.system_parms.number_of_projections_per_angle;
            xtal_pairs = zeros(2, num_of_angles * num_of_projs_per_angle);
            k = 1;
            for i=1:num_of_angles
                for j=1:num_of_projs_per_angle
                    pp = xp_first_angle(:, fix(nn / 2) - fix(num_of_projs_per_angle / 2) + j);
                    p0 = mod(pp(1) + i-1, xtal_num_trans_total);
                    p1 = mod(pp(2) + i-1, xtal_num_trans_total);
                    xtal_pairs(:,k) = [p0; p1];
                    k = k + 1;
                end
            end
            xtal_pairs = xtal_pairs + 1;
        end
        
        
               
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    methods(Access = 'public')
        
        
        
        function rp = getDefaultSinogramPlaneOrder(obj)
            %get default plane arrangement in terms of ring pair (different from Micxx-gram)
            %function rp = getDefaultSinogramPlaneOrder
            %
            nring = obj.getNumberOfCrystalRings();
            rp = zeros(2, nring*nring);
            offset = 0;
            for n = 1 : nring
                if n==1
                    rp(:,1:nring) = [1:nring; 1:nring];
                    offset = offset + nring;
                else
                    r_odd = [1:(nring-n+1); n:(nring)];
                    r_even = nring-r_odd+1;
                    nr = (nring-n+1)*2;
                    rp(:, offset + (1:2:nr)) = r_odd;
                    rp(:, offset + (2:2:nr)) = r_even;
                    offset = offset + nr;
                end
            end
        end
    end
    
    
    
    
    
    
    
        
    
    methods(Access = 'public')
        
        
        
        function prjs = doListModeForwardProjectionNonTOF(obj, image, image_size, voxel_size, lmdata)
            %do listmode-based forward projection for non-TOF case
            %function prjs = doListModeForwardProjectionNonTOF(image, image_size, voxel_size, lmdata)
            %
            %image = reshape(image, image_size);
            %
            if (~isa(lmdata, 'int16')) && (~isa(lmdata, 'uint16'))
                error('invalid data type: lmdata must be int16 or uint16!');
            end
            %
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            switch obj.system_parms.projector
                case 'Siddon'
                    tic
                    prjs = fproj_mt(image, image_size, voxel_size, tc, to, lmdata);
                    toc
                otherwise
                    error('invalid projector!');
            end
        end
        






        function prjs = doListModeForwardProjectionTOF(obj, image, image_size, voxel_size, lmdata)
            %do listmode-based forward projection for TOF case
            %function prjs = doListModeForwardProjectionNonTOF(image, image_size, voxel_size, lmdata)
            %
            image = reshape(image, image_size);
            %
            if ~isa(lmdata, 'int16')
                error('invalid data type: lmdata must be int16');
            end
            %
            if isempty(obj.system_parms.tof_info)
                error('unable to run TOF projector on non-TOF scanner!');
            end
            
            %
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            switch obj.system_parms.projector
                case 'Siddon'
                    tic
                    prjs = fproj_tof_mt(image, image_size, voxel_size, tc, to, lmdata, ...
                        obj.system_parms.tof_info);
                    toc
                otherwise
                    error('invalid projector!');
            end
        end
        


        
        function image = doListModeBackProjectionNonTOF(obj, prjs, image_size, voxel_size, lmdata)
            %do listmode-based backprojection for non-TOF case
            %function image = doListModeBackProjectionNonTOF(projs, image_size, voxel_size, lmdata)
            %
            %
            if (~isa(lmdata, 'int16')) && (~isa(lmdata, 'uint16'))
                error('invalid data type: lmdata must be int16 or uint16!');
            end
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            switch obj.system_parms.projector
                case 'Siddon'
                    tic
                    image = bproj_mt(prjs, image_size, voxel_size, tc, to, lmdata);
                    toc
                otherwise
                    error('invalid projector!');
            end
            image = reshape(image, image_size);
        end


     


        
        function image = doListModeBackProjectionTOF(obj, prjs, image_size, voxel_size, lmdata)
            %do listmode-based backprojection for TOF case
            %function image = doListModeBackProjectionTOF(projs, image_size, voxel_size, lmdata)
            %
            %
            if ~isa(lmdata, 'int16')
                error('invalid data type: lmdata must be int16 or uint16!');
            end
            
            if isempty(obj.system_parms.tof_info)
                error('unable to run TOF projector on non-TOF scanner!');
            end
            
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            switch obj.system_parms.projector
                case 'Siddon'
                    image = bproj_tof_mt(prjs, image_size, voxel_size, tc, to, lmdata, ...
                        obj.system_parms.tof_info);
                otherwise
                    error('invalid projector');
            end
            image = reshape(image, image_size);
            
        end
        
        
        
        





        function senimg = cal_senimg_uih(obj, plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, blockringdiff, num_blockrings, num_bins_axial_reduced, first_bed_ring, last_bed_ring)  
        
        	cryseff_wgap_inv = 1 ./ cryseff_wgap; 
        	plaeff_wgap_inv = 1 ./ plaeff_wgap; 
        	cryseff_wgap_inv(isinf(cryseff_wgap_inv)) = 0;
        	plaeff_wgap_inv(isinf(plaeff_wgap_inv)) = 0; 
        	
        	cryseff_wgap = cryseff_wgap_inv; 
        	plaeff_wgap = plaeff_wgap_inv; 
            
            
            %att_image = reshape(att_image, att_image_size);
            
            xp = int16(obj.getDefaultSinogramCrystalPairs());
            nring = obj.getNumberOfCrystalRings()
            
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            
            
            sino_xpairs = int16(obj.getDefaultSinogramCrystalPairs);
            size(sino_xpairs)
            
            
            nrad = obj.system_parms.number_of_projections_per_angle
            nang = obj.getDefaultNumberOfAngles
            
            num_bins_radial = nrad;
            num_bins_angular = nang;
            
            num_bins_sino = num_bins_radial * num_bins_angular
            
            crystal_array = obj.system_parms.crystal_array_size(1)
            num_blocks = obj.system_parms.number_of_detector_modules_transaxial
            
            num_crystals = crystal_array * num_blocks
            
            nocrypairs = zeros(num_crystals, num_crystals, 'uint32');
            
            
            
            for m=1:num_bins_sino
                
                nv = sino_xpairs(1, m);
                nu = sino_xpairs(2, m);
                
                nocrypairs(nv, nu) = m;
                nocrypairs(nu, nv) = m;
                
            end
            
            
            num_gaps = num_blockrings - 1
            
            nring_wogap = nring - num_gaps

            % [plane2] = my_ring_config_uih(nring_wogap);
            % no_plane2_pairs = zeros(nring, nring);
            % for nr = 1:nring_wogap^2
            %     nr1 = plane2(nr, 1);
            %     nr2 = plane2(nr, 2);
            %     noblock1 = floor((nr1-1)/num_bins_axial_reduced);
            %     noblock2 = floor((nr2-1)/num_bins_axial_reduced);
            %     no_plane2_pairs(nr1+noblock1, nr2+noblock2) = nr;
            % end
            % whos no_plane2_pairs
            % max(no_plane2_pairs(:))
            % min(no_plane2_pairs(:))



            
%              fprintf('calculating pure geometrical sensitivity for maximum block ring difference = %d...\n',  max_blockringdiff);           
            fprintf('calculating pure geometrical sensitivity for block ring difference = %d...\n',  blockringdiff);
            fprintf('image size: %d x %d x %d, voxel size: %f x %f x %f (mm^3)\n',  image_size, voxel_size);            
            %   M = obj.system_parms.crystal_array_size(1) * obj.system_parms.number_of_detector_modules_transaxial;

            
            i = int16(0);
            
            %    m = repmat(int16(0:nring-1),size(xp,2),1);
            
            
            m0=[];
            for noblock=0:num_gaps
                m0 = [m0, repmat(int16(noblock + noblock*num_bins_axial_reduced : noblock-1 + (noblock+1)*num_bins_axial_reduced ), 1, 1)];
            end
            
            
            m=[];
            for noblock=0:num_gaps
                m = [m, repmat(int16(noblock + noblock*num_bins_axial_reduced : noblock-1 + (noblock+1)*num_bins_axial_reduced ), size(xp,2), 1)];
            end
            
            
            % whos 
           
            lmdata = repmat(i, 5, size(xp,2) * nring_wogap);
            lmdata(1,:) = repmat(xp(1,:)-1, 1, nring_wogap);
            lmdata(3,:) = repmat(xp(2,:)-1, 1, nring_wogap);
            lmdata(4,:) = m(:);
            
            senimg = 0;
            
            num_xtals_wgap = num_bins_axial_reduced + 1;
            
            
            num_bins_angular_reduced = crystal_array;
            
            
            %for n = 1: nring_wogap %  589  %      1 : nring_wogap
            for n = first_bed_ring : last_bed_ring %  589  %      1 : nring_wogap
                
                fprintf('processing ring without gap   #%d ...\n', n);


                if (  ( ( ceil(nring - m0(n))/num_xtals_wgap ) >= blockringdiff) ||  ( ( ceil(m0(n))/num_xtals_wgap ) >= blockringdiff ) )
                
                    lmdata(2,:) = int16(m0(n));
                    
                    noblockring1 = floor(single(lmdata(2,:))/num_xtals_wgap);
                    noblockring2 = floor(single(lmdata(4,:))/num_xtals_wgap);

    %                 max_noblockring1 = max(noblockring1)
    %                 min_noblockring1 = min(noblockring1)
    %                 max_noblockring2 = max(noblockring2)
    %                 min_noblockring2 = min(noblockring2)

                    noblcokringdifference = uint32(abs(noblockring1 - noblockring2));
                    
    %                 max_noblcokringdifference = max(noblcokringdifference)
    %                 min_noblcokringdifference = min(noblcokringdifference)

                    clear  noblockring1  noblockring2
                                    


                    % idx1 = int16(noblcokringdifference == blockringdiff);
                    idx2 = find(noblcokringdifference == blockringdiff);

                    % whos lmdata idx1 idx2

                    lm0 = lmdata(:, idx2);  


                    % lm0_noblockring1 = floor(single(lm0(2,:))/num_xtals_wgap);
                    % lm0_noblockring2 = floor(single(lm0(4,:))/num_xtals_wgap); 
                    % lm0 = lmdata(:, lmdata(2,:)<=lmdata(4,:) & int16(noblcokringdifferenceM == blockringdiff) & int16(noblockring1M>=0) & int16(noblockring2M>=0) & int16(noblockring1M<=6) & int16(noblockring2M<=6));   

                    % whos no_plane2_pairs cryseff_wgap lm0 lmdata idx1 idx2 lm0_noblockring*

                    % max(lm0(4,:))
                    % min(lm0(4,:))
                    % max(lm0(3,:))
                    % min(lm0(3,:))

                    % max(lm0(2,:)+1)
                    % min(lm0(2,:)+1)
                    % max(lm0(4,:)+1)
                    % min(lm0(4,:)+1)



                    % linearInd_1 = sub2ind(size(cryseff_wgap), uint32(lm0(2,:)+1), uint32(lm0(1,:)+1) );
                    % lm_cryseff_wgap_1 = cryseff_wgap(linearInd_1);
                    % linearInd_2 = sub2ind(size(cryseff_wgap), uint32(lm0(4,:)+1), uint32(lm0(3,:)+1) );
                    % lm_cryseff_wgap_2 = cryseff_wgap(linearInd_2);

                    % linearInd_3 = sub2ind(size(no_plane2_pairs), uint32(lm0(2,:)+1), uint32(lm0(4,:)+1) );
                    % % max(linearInd_3)
                    % % min(linearInd_3)
                    % lm_plaeff = plaeff(no_plane2_pairs(linearInd_3)); 

                    % % whos lm_cryseff_wgap_1 lm_cryseff_wgap_2 lm_plaeff

                    linearInd_1 = sub2ind(size(cryseff_wgap), uint32(lm0(2,:)+1), uint32(lm0(1,:)+1) );
                    lm_cryseff_wgap_1 = cryseff_wgap(linearInd_1);
                    linearInd_2 = sub2ind(size(cryseff_wgap), uint32(lm0(4,:)+1), uint32(lm0(3,:)+1) );
                    lm_cryseff_wgap_2 = cryseff_wgap(linearInd_2);
                    linearInd_3 = sub2ind(size(plaeff_wgap), uint32(lm0(2,:)+1), uint32(lm0(4,:)+1) );
                    lm_plaeff_wgap = plaeff_wgap(linearInd_3); 


                    % v1: seems wrong plane efficiency: should be inversed
                    lm_nrm = double(lm_cryseff_wgap_1) .* double(lm_cryseff_wgap_2) .* double(lm_plaeff_wgap);
                    % v2: e1*e2/planeff
                    % lm_nrm = double(lm_cryseff_wgap_1) .* double(lm_cryseff_wgap_2) ./ double(lm_plaeff_wgap);
                    % v3: 1/(e1*e2*planeff)
                    %lm_nrm = 1 ./ (double(lm_cryseff_wgap_1) .* double(lm_cryseff_wgap_2) .* double(lm_plaeff_wgap));
            
                    % raw_data = double(lm_cryseff_wgap_1) .* double(lm_cryseff_wgap_2) .* double(lm_plaeff');

                    % whos raw_data

                    % lm0(:, 1:10)
                    % raw_data(1:10)

                    clear  noblockring1  noblockring2   noblcokringdifference  ringdifference ...
                          linearInd   linearInd2 index_angular_bin   index_angular_rebin ...
                          index_radial_no   index_angular_no    index_axial_xp_no   index_blrgdiff_no
                                      


                    % t0 = size(lm_nrm)
                    
                    %att_image = 0 .* att_image; 
                    prjs_att = fproj_mt(att_image, att_image_size, att_voxel_size, tc, to, lm0);
                    
                    t = size(prjs_att);

                    
                    ratio_MC_Ana = reshape(lm_nrm, t(1), t(2));
                    

                    clear lm_nrm
                    
                    s0 = obj.doListModeBackProjectionNonTOF( double(ratio_MC_Ana) .* exp(-prjs_att), image_size, voxel_size, lm0);

                    
                    senimg = senimg + s0;
                    
                    
                    clear lm_nrm  prjs_att  ratio_MC_Ana  lm0 s0
                                  
                    
                end

                 
            end

                                
            whos

            
            senimg = reshape(senimg, image_size);
                
             

        end%function cal_senimg_uih  
        
        function senimg = cal_senimg_single_bed(obj, plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, num_blockrings, num_bins_axial_reduced, bedStartRing, bedEndRing)
            %   bedStartRing: ring number of the start of the bed position (from 1 to 672)
            %   bedEndRing: ring number of end  of bed
            senimg = 0;
            
            bedStartRing
            bedEndRing

            % calculated inverted crystal efficiencies and plane efficiencies
            cryseff_wgap_inv = 1 ./ cryseff_wgap;
            plaeff_wgap_inv = 1 ./ plaeff_wgap; 
            cryseff_wgap_inv(isinf(cryseff_wgap_inv)) = 0;
            plaeff_wgap_inv(isinf(plaeff_wgap_inv)) = 0; 
            
            cryseff_wgap = cryseff_wgap_inv; 
            plaeff_wgap = plaeff_wgap_inv; 
            
            nring = obj.getNumberOfCrystalRings()  %679, includes gaps
            
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            
            
            sino_xpairs = int16(obj.getDefaultSinogramCrystalPairs);    %2x230580 matrix (rows x columns)
            %  size(sino_xpairs);

            
            nrad = obj.system_parms.number_of_projections_per_angle;
            nang = obj.getDefaultNumberOfAngles;
            
            num_bins_radial = nrad     %549
            num_bins_angular = nang    %420
            
            num_bins_sino = num_bins_radial * num_bins_angular                     % 230580
            
            crystal_array = obj.system_parms.crystal_array_size(1)                 % 35
            num_blocks = obj.system_parms.number_of_detector_modules_transaxial    % 24
            
            num_crystals = crystal_array * num_blocks                              % 840
            
            nocrypairs = zeros(num_crystals, num_crystals, 'uint32');         

            
            num_gaps = num_blockrings - 1                                          % 7
            
            nring_wogap = nring - num_gaps                                     % 672


            
%              fprintf('calculating pure geometrical sensitivity for maximum block ring difference = %d...\n',  max_blockringdiff);           
         %   fprintf('calculating pure geometrical sensitivity for block ring difference = %d...\n',  blockringdiff);
            fprintf('image size: %d x %d x %d, voxel size: %f x %f x %f (mm^3)\n',  image_size, voxel_size);            
            %   M = obj.system_parms.crystal_array_size(1) * obj.system_parms.number_of_detector_modules_transaxial;

            
            i = int16(0);
                        
            
            m0=[];

            % creates an array m0 with numbers from 0 to 678 but omitting the axial gaps, i.e. 84, 169, 254, 339, ...
            % int16(a : b) creates a list of integers (row vector) starting with a incrementing in steps of 1 and finishing with b as last entry
            for noblock=0:num_gaps
                m0 = [m0, repmat(int16(noblock + noblock*num_bins_axial_reduced : noblock-1 + (noblock+1)*num_bins_axial_reduced ), 1, 1)]; %repmat() has no effect in this case
            end
            
            
            
            m=[];
            % creates an array with entries from 0 to 678 omitting gaps (like above) and does that for every sinogram bin (230580)
            % -> rows: Numbers from 0 to 678 omitting gaps
            % -> 230580 rows in total
            for noblock=0:num_gaps
                m = [m, repmat(int16(noblock + noblock*num_bins_axial_reduced : noblock-1 + (noblock+1)*num_bins_axial_reduced ), size(sino_xpairs,2), 1)];
            end
            %whos
            
           
            lmdata = repmat(i, 5, size(sino_xpairs,2) * nring_wogap);   % array filled with value i at each position; has 5 rows, each with 230580 * 672 entries
            lmdata(1,:) = repmat(sino_xpairs(1,:)-1, 1, nring_wogap);   % decrement all elements in 1st row of sino_xpairs by 1, write row nring_wogap times in a row -> save in lmdata
            lmdata(3,:) = repmat(sino_xpairs(2,:)-1, 1, nring_wogap);   % decrement all elements in second row of sino_xpairs by 1, write row nring_wogap times in a row-> save 
            lmdata(4,:) = m(:); % 4th row of lmdata filled with elements in m, row by row, i.e. 230580 elements per row, 672 times:
                                % 0, 0, 0, ..., 1, 1, 1, ..., 679, 679, 679, .... 679; each number appears 230580 times

            
            num_xtals_wgap = num_bins_axial_reduced + 1;    % num_xtals_wgap = 84 +1 = 85;
            
            
            num_bins_angular_reduced = crystal_array;       % 35
            
            for n = bedStartRing : bedEndRing % nring_wogap 
                
                fprintf('processing ring without gap   #%d ...\n', n);


                %if (  ( ( ceil(nring - m0(n))/num_xtals_wgap ) >= blockringdiff) ||  ( ( ceil(m0(n))/num_xtals_wgap ) >= blockringdiff ) )
                
                    lmdata(2,:) = int16(m0(n)); % fill 2nd row of lmdata with current crystal ring: m0(1) = 0, m0(2) = 1, etc
                %{    
                    noblockring1 = floor(single(lmdata(2,:))/num_xtals_wgap);   %fill array noblockring1 with current unit from 0 to 7
                    noblockring2 = floor(single(lmdata(4,:))/num_xtals_wgap);   %fill array noblockring2 with unit of all rings


                    noblcokringdifference = uint32(abs(noblockring1 - noblockring2));
                        % calculates the unit difference between current ring and any other ring

                    clear  noblockring1  noblockring2
                                    
                    idx2 = find(noblcokringdifference == blockringdiff);    % get ALL positions in array noblcokringdifference with correct unit difference (i.e. equals blockringdiff) 
                                                                            % and write them to array indx2

                %}
                    lm0 = [];
                    for r = bedStartRing : bedEndRing
                        a = (r-1)*size(sino_xpairs,2)+1;
                        b = (r)*size(sino_xpairs,2);
                        lm0 = [lm0, lmdata(:,[a:b])];
                    end
                    %lm0 = lmdata(:, idx2);  %saves every column of lmdata with matching unit difference and saves it to lm0.


                    linearInd_1 = sub2ind(size(cryseff_wgap), uint32(lm0(2,:)+1), uint32(lm0(1,:)+1) );
                    lm_cryseff_wgap_1 = cryseff_wgap(linearInd_1);
                    linearInd_2 = sub2ind(size(cryseff_wgap), uint32(lm0(4,:)+1), uint32(lm0(3,:)+1) );
                    lm_cryseff_wgap_2 = cryseff_wgap(linearInd_2);
                    linearInd_3 = sub2ind(size(plaeff_wgap), uint32(lm0(2,:)+1), uint32(lm0(4,:)+1) );
                    lm_plaeff_wgap = plaeff_wgap(linearInd_3); 


                    lm_nrm = double(lm_cryseff_wgap_1) .* double(lm_cryseff_wgap_2) .* double(lm_plaeff_wgap);


                    clear  noblockring1  noblockring2   noblcokringdifference  ringdifference ...
                          linearInd   linearInd2 index_angular_bin   index_angular_rebin ...
                          index_radial_no   index_angular_no    index_axial_xp_no   index_blrgdiff_no
                                      

                    prjs_att = fproj_mt(att_image, att_image_size, att_voxel_size, tc, to, lm0);
                    
                    t = size(prjs_att);

                    
                    ratio_MC_Ana = reshape(lm_nrm, t(1), t(2));
                    

                    clear lm_nrm
                    
                    s0 = obj.doListModeBackProjectionNonTOF( double(ratio_MC_Ana) .* exp(-prjs_att), image_size, voxel_size, lm0);

                    
                    senimg = senimg + s0;
                    
                    
                    clear lm_nrm  prjs_att  ratio_MC_Ana  lm0 s0
                                  
                    
                %end

        %}   
            end% for loop

            senimg = reshape(senimg, image_size);




            


            %whos
        end


        function senimg = cal_senimg_multi_bed(obj, plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, num_blockrings, num_bins_axial_reduced, bedStartRing, numRingsPerBed, numBeds, overlap)
            %   bedStartRing: ring number of the start of the FIRST bed position (from 1 to 672)
            %   numRingsPerBed:     number of rings in each bed positions
            %   numBeds: number of overlapping bed positions
            %   overlap: fraction of rings that overlap for two adjacent bed positions

            % calculate start and end ring for each bed positions
            StartRings = zeros(numBeds, 1); % start ring of each bed position
            StartRings(1,1) = bedStartRing;
            EndRings = zeros(numBeds, 1);   % end ring of each bed position
            EndRings(1,1) = bedStartRing+numRingsPerBed-1;
            overlappingRings = round(numRingsPerBed*overlap);  % number of rings that are overlapping
            for n = 2 : numBeds
                StartRings(n,1) = StartRings(n-1,1) + numRingsPerBed - overlappingRings;
                EndRings(n,1) = StartRings(n,1) + numRingsPerBed-1;
            end
            
            senimg = 0;

            % calculated inverted crystal efficiencies and plane efficiencies
            cryseff_wgap_inv = 1 ./ cryseff_wgap;
            plaeff_wgap_inv = 1 ./ plaeff_wgap; 
            cryseff_wgap_inv(isinf(cryseff_wgap_inv)) = 0;
            plaeff_wgap_inv(isinf(plaeff_wgap_inv)) = 0; 
            
            cryseff_wgap = cryseff_wgap_inv; 
            plaeff_wgap = plaeff_wgap_inv; 
            
            nring = obj.getNumberOfCrystalRings();  %679, includes gaps
            
            tc = obj.getCrystalTransaxialLocations();
            to = obj.getCrystalRingOffsets();
            
            
            sino_xpairs = int16(obj.getDefaultSinogramCrystalPairs);    %2x230580 matrix (rows x columns)
            %  size(sino_xpairs);
            
            nrad = obj.system_parms.number_of_projections_per_angle;
            nang = obj.getDefaultNumberOfAngles;
            
            num_bins_radial = nrad;     %549
            num_bins_angular = nang;    %420
            
            num_bins_sino = num_bins_radial * num_bins_angular;                     % 230580
            
            crystal_array = obj.system_parms.crystal_array_size(1);                 % 35
            num_blocks = obj.system_parms.number_of_detector_modules_transaxial;    % 24
            
            num_crystals = crystal_array * num_blocks;                              % 840
            
            nocrypairs = zeros(num_crystals, num_crystals, 'uint32');         

            
            num_gaps = num_blockrings - 1;                                          % 7
            
            nring_wogap = nring - num_gaps;                                         % 672


            
%              fprintf('calculating pure geometrical sensitivity for maximum block ring difference = %d...\n',  max_blockringdiff);           
         %   fprintf('calculating pure geometrical sensitivity for block ring difference = %d...\n',  blockringdiff);
            fprintf('image size: %d x %d x %d, voxel size: %f x %f x %f (mm^3)\n',  image_size, voxel_size);            
            %   M = obj.system_parms.crystal_array_size(1) * obj.system_parms.number_of_detector_modules_transaxial;

            
            i = int16(0);
                        
            
            m0=[];

            % creates an array m0 with numbers from 0 to 678 but omitting the axial gaps, i.e. 84, 169, 254, 339, ...
            % int16(a : b) creates a list of integers (row vector) starting with a incrementing in steps of 1 and finishing with b as last entry
            for noblock=0:num_gaps
                m0 = [m0, repmat(int16(noblock + noblock*num_bins_axial_reduced : noblock-1 + (noblock+1)*num_bins_axial_reduced ), 1, 1)]; %repmat() has no effect in this case
            end
            
            
            
            m=[];
            % creates an array with entries from 0 to 678 omitting gaps (like above) and does that for every sinogram bin (230580)
            % -> rows: Numbers from 0 to 678 omitting gaps
            % -> 230580 rows in total
            for noblock=0:num_gaps
                m = [m, repmat(int16(noblock + noblock*num_bins_axial_reduced : noblock-1 + (noblock+1)*num_bins_axial_reduced ), size(sino_xpairs,2), 1)];
            end
            %whos
            
           
            lmdata = repmat(i, 5, size(sino_xpairs,2) * nring_wogap);   % array filled with value i at each position; has 5 rows, each with 230580 * 672 entries
            lmdata(1,:) = repmat(sino_xpairs(1,:)-1, 1, nring_wogap);   % decrement all elements in 1st row of sino_xpairs by 1, write row nring_wogap times in a row -> save in lmdata
            lmdata(3,:) = repmat(sino_xpairs(2,:)-1, 1, nring_wogap);   % decrement all elements in second row of sino_xpairs by 1, write row nring_wogap times in a row-> save 
            lmdata(4,:) = m(:); % 4th row of lmdata filled with elements in m, row by row, i.e. 230580 elements per row, 672 times:
                                % 0, 0, 0, ..., 1, 1, 1, ..., 679, 679, 679, .... 679; each number appears 230580 times

            
            num_xtals_wgap = num_bins_axial_reduced + 1;    % num_xtals_wgap = 84 +1 = 85;
            
            
            num_bins_angular_reduced = crystal_array;       % 35
            
            latestEndRing = StartRings(1,1);    %latestEndRing saves the ring number of the last bed position

            for nbed = 1 : numBeds % loop over all bed positions
                fprintf('processing bed position   #%d ...\n', nbed);

                for n = latestEndRing : EndRings(nbed,1) % loop over all rings of current bed position that haven't been treated in the last position (i.e. those rings exceeding the overlap) 
                    
                    fprintf('processing ring without gap   #%d ...\n', n);
                    
                        lmdata(2,:) = int16(m0(n)); % fill 2nd row of lmdata with current crystal ring: m0(1) = 0, m0(2) = 1, etc

                        lm0 = [];
                        for r = StartRings(nbed,1) : EndRings(nbed,1)
                            a = (r-1)*size(sino_xpairs,2)+1;
                            b = (r)*size(sino_xpairs,2);
                            lm0 = [lm0, lmdata(:,[a:b])];
                        end

                        linearInd_1 = sub2ind(size(cryseff_wgap), uint32(lm0(2,:)+1), uint32(lm0(1,:)+1) );
                        lm_cryseff_wgap_1 = cryseff_wgap(linearInd_1);
                        linearInd_2 = sub2ind(size(cryseff_wgap), uint32(lm0(4,:)+1), uint32(lm0(3,:)+1) );
                        lm_cryseff_wgap_2 = cryseff_wgap(linearInd_2);
                        linearInd_3 = sub2ind(size(plaeff_wgap), uint32(lm0(2,:)+1), uint32(lm0(4,:)+1) );
                        lm_plaeff_wgap = plaeff_wgap(linearInd_3); 


                        lm_nrm = double(lm_cryseff_wgap_1) .* double(lm_cryseff_wgap_2) .* double(lm_plaeff_wgap);


                        clear  noblockring1  noblockring2   noblcokringdifference  ringdifference ...
                              linearInd   linearInd2 index_angular_bin   index_angular_rebin ...
                              index_radial_no   index_angular_no    index_axial_xp_no   index_blrgdiff_no
                                          

                        prjs_att = fproj_mt(att_image, att_image_size, att_voxel_size, tc, to, lm0);
                        
                        t = size(prjs_att);

                        
                        ratio_MC_Ana = reshape(lm_nrm, t(1), t(2));
                        

                        clear lm_nrm
                        
                        s0 = obj.doListModeBackProjectionNonTOF( double(ratio_MC_Ana) .* exp(-prjs_att), image_size, voxel_size, lm0);

                        
                        senimg = senimg + s0;
                        
                        
                        clear lm_nrm  prjs_att  ratio_MC_Ana  lm0 s0
                                      
                        
                    %end

            %}   
                end% for loop over rings

                latestEndRing = EndRings(nbed,1);

            end%for loop over beds

            senimg = reshape(senimg, image_size);




            


            %whos
        end

    end%methods
    
    
    
    
    
    
    
    
    methods
        
        
        
        function h = l2h(obj, lm)
            %convert listmode data to histograms
            %function h = l2h(lmdata)
            if (size(lm,1)~= 4 && size(lm,1)~= 5)
                error('unknown listmode data format, should be 4xN or 5xN (TOF) (currently not support DOI bins)');
            end
            
            M = obj.system_parms.crystal_array_size(1) * obj.system_parms.number_of_detector_modules_transaxial;
            nring = obj.getNumberOfCrystalRings;
            
            % remove unexpected prompts
            ii=((lm(1,:)<M) & (lm(1,:)>=0)) & ...
                ((lm(3,:)<M) & (lm(3,:)>=0)) & ...
                ((lm(2,:)<nring) & (lm(2,:)>=0)) & ...
                ((lm(4,:)<nring) & (lm(4,:)>=0));
            
            if length(ii) ~= size(lm, 2)
                warning('invalid listmode data found, # = %d', size(lm,2)-length(ii));
            end
            
            lm = lm(:,ii);
            idx = sub2ind([M, M, nring, nring], ...
                double(lm(1,:))+1, double(lm(3,:))+1, double(lm(2,:))+1, double(lm(4,:))+1);
            
            h_M = M
            h_nring = nring

            h = zeros(M, M, nring, nring);
            ss = sort(idx, 'ascend');

            % idx(1:10)
            % idx(end-10:end)

            % whos idx
            clear('idx');

            % ss(1:10)
            % ss(end-10:end)

            nz=0;
            c=1;


            while (~isempty(ss))
                
                dd = diff(ss);

                % dd(1:10)
                % dd(end-10:end)

                iii = [1, find(dd>0)+1]; 

                % iii(1:10)
                % iii(end-10:end)

                nz=nz + length(iii);

                % length(iii)
                % nz

                h(ss(iii)) = h(ss(iii)) + 1;
                
                ss(iii)=-1;

                ss=ss(ss>0);

                % % if length(ss) > 10
                %     ss(1:10)
                %     ss(end-10:end)
                % % else
                % %     ss(1:length(ss))
                % % end

                if 0
                    fprintf('#%d: %d (%d)\n', c, length(iii), size(lm,2));
                end

                c=c+1;

                c
                
            end
            fprintf('nz=%d\n', nz);
        end
        
        
        
        

        function s = h2s(obj, h, flag)
        %convert histogram to singoram (h: MxMxnringxnring)
        %function s = h2s(h)
            %fprintf('%d\n', nargin);
            if nargin < 3
                flag = false;
            end
            xp = obj.getDefaultSinogramCrystalPairs;
            idx = sub2ind([size(h,1), size(h,2)], xp(1,:), xp(2,:));
            s = zeros(size(xp,2), size(h,3), size(h,4));
            for n = 1 : size(h,4)
                for m = 1 : size(h,3)
                    if flag % alread joint!
                        hh = h(:,:,m,n);
                    else
                        hh = h(:,:,m,n) + h(:,:,n,m)';
                    end
                    s(:, m, n) = hh(idx);
                end
            end
            s = reshape(s, obj.system_parms.number_of_projections_per_angle, ...
                        obj.getDefaultNumberOfAngles, size(h,3), size(h,4));
            disp('NOTE: natural order!');
        end




        
        
        
    end
    
    
end
