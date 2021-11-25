function plane_eff = compute_plaeff_map(lm_file, attn_coeff_file, cryseff_file)

% script computes plane efficiencies based on homogenous cylinder data
% also computes mean attenuation factors for every plane

%% fixed parameters
num_crys_ax = 679;      % gaps included
num_crys_trans = 840;

plane_eff = zeros(num_crys_ax, num_crys_ax);
entries = zeros(num_crys_ax, num_crys_ax);
entries_crystal_map = zeros(num_crys_ax, num_crys_trans);
mean_attn_per_plane = zeros(num_crys_ax, num_crys_ax);
mean_attn_per_crystal = zeros(num_crys_ax, num_crys_trans);

fid_crys = fopen(cryseff_file, 'rb');
crys_eff = fread(fid_crys, inf, 'float');
crys_eff = reshape(crys_eff, num_crys_ax, num_crys_trans);


%% open files in 10 portions
% necessary to prevent matlab from crashing when reading in too large files
% at once; rbayerlein; 10/28/2021
for p = 1:10
    fprintf('now processing %d\n', p);
    fid = fopen(lm_file, 'rb');
    s = dir(lm_file);
    num_events = double(s.bytes)/10;
    num_to_read = floor(num_events/10);

    fseek(fid, num_to_read*10*(p-1), 'bof');    % skip (p-1) times number of events with 10 byte each
    lm_data = fread(fid, [5 num_to_read], '*int16');    % read in 5 values per column of size int16 (=2bytes) each

    lm_data = reshape(lm_data, 5, num_to_read);

    fid_attn = fopen(attn_coeff_file, 'rb');
    fseek(fid_attn, num_to_read*4*(p-1), 'bof');% skip (p-1) times number of events with 4 byte each
    attn_data = fread(fid_attn, num_to_read, 'float');

    %% main program

    for evt = 1 : num_to_read
        plane_eff(lm_data(2,evt)+1,lm_data(4,evt)+1) = plane_eff(lm_data(2,evt)+1,lm_data(4,evt)+1) + 1*crys_eff(lm_data(2,evt)+1,lm_data(1,evt)+1)*crys_eff(lm_data(4,evt)+1,lm_data(3,evt)+1)/attn_data(evt);
        entries(lm_data(2,evt)+1,lm_data(4,evt)+1) = entries(lm_data(2,evt)+1,lm_data(4,evt)+1) + 1;
        entries_crystal_map(lm_data(2,evt)+1,lm_data(1,evt)+1) = entries_crystal_map(lm_data(2,evt)+1,lm_data(1,evt)+1) + 1;
        entries_crystal_map(lm_data(4,evt)+1,lm_data(3,evt)+1) = entries_crystal_map(lm_data(4,evt)+1,lm_data(3,evt)+1) + 1;
        mean_attn_per_plane(lm_data(2,evt)+1,lm_data(4,evt)+1) = mean_attn_per_plane(lm_data(2,evt)+1,lm_data(4,evt)+1)+attn_data(evt);
        mean_attn_per_crystal(lm_data(2,evt)+1,lm_data(1,evt)+1) = mean_attn_per_crystal(lm_data(2,evt)+1,lm_data(1,evt)+1) + attn_data(evt);
        mean_attn_per_crystal(lm_data(4,evt)+1,lm_data(3,evt)+1) = mean_attn_per_crystal(lm_data(4,evt)+1,lm_data(3,evt)+1) + attn_data(evt);
    end
    fclose(fid);
    fclose(fid_attn);

end
mean_attn_per_plane = mean_attn_per_plane ./ entries;
mean_attn_per_crystal = mean_attn_per_crystal ./ entries_crystal_map;

section = mean_attn_per_crystal([640 650],:);
mean(mean(section))
section_light = mean_attn_per_crystal([620 630],:);
mean(mean(section_light))

% mean(mean(entries(entries>0)))

%imshow(plane_eff, []);
%colorbar
figure;
imshow(mean_attn_per_plane, []);
colorbar;
figure;
imshow(entries_crystal_map, []);
colorbar;
figure;
imshow(mean_attn_per_crystal, []);
colorbar;