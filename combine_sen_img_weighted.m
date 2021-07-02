function combine_sen_img(senimg_name, num_beds, rings_per_bed, bedStartRing, overlap_percent)

dim_x = 239;
dim_y = dim_x;
dim_z = 679;

sen_img = zeros (dim_x, dim_y, dim_z);
sen_img = sen_img(:);

overlap = overlap_percent/100

lastBedEndRing = bedStartRing + (1+ (num_beds-1)*(1-overlap))*rings_per_bed-1

slice_weighting = zeros(dim_z,1);

for r = 1 : num_beds
    start = bedStartRing + (r-1)*(1-overlap)*rings_per_bed
    finish = start + rings_per_bed-1
    for s = start : finish
        slice_weighting(s,1) = slice_weighting(s,1)+1;
    end
    senimg_name_temp = senimg_name(1:strfind(senimg_name, '.sen_img')); 
    senimg_name_temp = [senimg_name_temp, num2str(start), '_', num2str(finish), '.sen_img'];
    ss = ['now reading in ', senimg_name_temp];
    disp(ss);

    fid = fopen(senimg_name_temp, 'rb');
    sen_temp = fread(fid, inf, 'float');
    sen_img = sen_img + sen_temp;
    fclose(fid);
end

sen_img = reshape(sen_img, dim_x, dim_y, dim_z);
for s = bedStartRing : lastBedEndRing
    if slice_weighting(s,1) ~= 0
        sen_img(:,:,s) = sen_img(:,:,s)/slice_weighting(s,1);
    end
end

out_name = senimg_name(1:strfind(senimg_name, '.sen_img'));
out_name = [out_name, 'combine_weighted_', num2str(num_beds), 'beds_', num2str(rings_per_bed), 'rings_from', num2str(bedStartRing), '_', num2str(overlap_percent), '%overlap', '.sen_img'];

ss=['output file name ' , out_name]
disp(ss);



slice = reshape(sen_img(:,round(239/2),:), [239, 679]);
imshow(slice, []);
colorbar;


%% fill rest of FOV with high values to get rid of reconstruction artifacts
%calculate lower end of field of view as ring number inlcuding gap
num_gaps_start = fix(bedStartRing/84)
start_ring = bedStartRing + num_gaps_start

%calculate upper end of FOV including gap
final_ring = bedStartRing + rings_per_bed + (1-overlap)*rings_per_bed*(num_beds-1)-1
num_gaps_end = fix(final_ring/84)
final_ring = final_ring + num_gaps_end
if rem(final_ring-num_gaps_end,84) == 0
    final_ring = final_ring-1
end

highVal = realmax('single')

img_MOD = zeros(dim_x, dim_y, dim_z);

for sl = 1:dim_z
    if sl < start_ring || sl > final_ring
        for x = 1 : dim_x   %
            for y = 1 : dim_y
                img_MOD(x,y,sl) = highVal;
            end
        end
    else
        for x = 1 : dim_x   %
            for y = 1 : dim_y
                img_MOD(x,y,sl) = sen_img(x,y,sl);
            end
        end
    end
end

fid_out = fopen(out_name, 'w');
fwrite(fid_out, img_MOD, 'float');
fclose(fid_out);


%fid_out_realName = fopen(senimg_name, 'w');
%fwrite(fid_out_realName, img_MOD, 'float');
%fclose(fid_out_realName);

end