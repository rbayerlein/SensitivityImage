function plane_eff = get_plaeff(nc_path, num_beds, first_bed_ring, rings_per_bed, overlap);

num_crys_ax = 672; 
last_ring = first_bed_ring+rings_per_bed+(num_beds-1)*(rings_per_bed-overlap) -1
if last_ring > num_crys_ax
    disp('Number of beds exceeds total length of EXPLORER field of view. Abort.')
    return;
end

if (overlap < 0 || overlap > rings_per_bed-1)
    disp('Bed overlap is not valid. Abort.')
    return;
end

%% open and read plane efficiency from file
fid_mich = fopen('michel_lut_672x672', 'rb'); 
if fid_mich < 0
	disp('Could not open michelogram LUT...'); 
	fid_mich = fopen('michel_lut_672x672'); 
	if fid_mich < 0
 		disp('Could not open michelogram LUTddd, quit');
		return;
	end
end

michel_lut_672x672 = fread(fid_mich, inf, 'int'); 
michel_lut_672x672 = reshape(michel_lut_672x672, num_crys_ax, num_crys_ax); 
fclose(fid_mich); 

plane_eff = ones(num_crys_ax,num_crys_ax); 

fid_nc = fopen(nc_path, 'rb'); 

if fid_nc < 0
	disp('Could not open .nc file, quit'); 
	return; 
end

fseek(fid_nc, 2079453*4, 'bof'); 

plane_eff_1 = fread(fid_nc, num_crys_ax*num_crys_ax, 'float');
fclose(fid_nc);


for k = 1:num_crys_ax
	for j = 1:num_crys_ax
		plane_eff(k,j) = plane_eff_1(michel_lut_672x672(k,j)+1); 
        plane_eff(k,j) = 1/plane_eff(k,j);  % inverse plane efficiency map
	end
end

%% calculate the average of the ORIGINAL inverted plane efficiency

average = 0;
count =0;
for x = 1 :num_crys_ax
    for y = 1:num_crys_ax
        if plane_eff(x,y) <10e30 && plane_eff(x,y) ~= 0 && plane_eff(x,y) ~=1
            average = average + plane_eff(x,y);
            count = count + 1;
        end
    end
end
average = average/count;
ss = ['average of ORIGINAL inverted plane efficiency map: ', num2str(average)];
disp(ss);

%% add plane efficiencies from individual bed positions
plane_eff_temp = zeros(num_crys_ax,num_crys_ax); 

for pos = 1 : num_beds
    startRing = first_bed_ring + (pos-1)*(rings_per_bed-overlap);
    endRing = startRing + rings_per_bed - 1;
    for ax = 1 : num_crys_ax
        if (ax >=startRing && ax <=endRing)
            for trans = 1 : num_crys_ax
                if (trans >=startRing && trans <=endRing)
                    plane_eff_temp(ax, trans) = plane_eff_temp(ax, trans) + plane_eff(ax, trans);
                end
            end
        end
    end            
end

%% calculate the average of the MODIFIED inverted plane efficiency

average_mod = 0;
count_mod =0;
for x = 1 :num_crys_ax
    for y = 1:num_crys_ax
        if plane_eff_temp(x,y) <10e30 && plane_eff_temp(x,y) ~= 0 && plane_eff_temp(x,y) ~=1
            average_mod = average_mod + plane_eff_temp(x,y);
            count_mod = count_mod + 1;
        end
    end
end
average_mod = average_mod/count_mod;
ss = ['average of MODIFIED inverted plane efficiency map: ', num2str(average_mod)];
disp(ss);

%% Calculate and apply scaling factor
scaling_factor = average/average_mod;
ss = ['scaling factor: ', num2str(scaling_factor)];
disp(ss);

for i = 1:num_crys_ax
    for j = 1:num_crys_ax
        if plane_eff_temp(i,j) < 10e30 && plane_eff_temp(i,j) ~=0 && plane_eff_temp(i,j) ~=1
            plane_eff_temp(i,j) = plane_eff_temp(i,j)*scaling_factor;
        end
    end
end

 %% invert the plane efficiencies back

 lowest_value = 0;
for ax = 1 : num_crys_ax
    for trans = 1 : num_crys_ax
        if plane_eff_temp(ax, trans) == 0
            plane_eff_temp(ax, trans) = lowest_value;
        end
        plane_eff(ax, trans) = 1/plane_eff_temp(ax, trans);
    end
end

% imshow(plane_eff, [0, 2])
% colorbar;

 %% add axial gap
%num_crys_ax_wgap = 679; 
%num_crys_ax_unit = 84;
 
%plane_eff_wgap = zeros(num_crys_ax_wgap, num_crys_ax_wgap); 
%for k = 0:(num_crys_ax-1)
%	for kk = 0:(num_crys_ax-1)
%		ax_new1 = floor(k/num_crys_ax_unit) + k; 
%		ax_new2 = floor(kk/num_crys_ax_unit) + kk; 
%		plane_eff_wgap(ax_new1+1, ax_new2+1) = plane_eff(k+1, kk+1); 
%	end
%end

%figure
%imagesc(plane_eff_wgap)

 