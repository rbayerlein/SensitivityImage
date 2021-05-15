function combine_sen_img(senimg_name, num_beds, ring_per_bed, bedStartRing)

sen_img = zeros (239,239,679);
sen_img = sen_img(:);

bed_shift = ring_per_bed/2

for r = 1 : num_beds
    start = bedStartRing + (r-1)*bed_shift
    finish = start + ring_per_bed-1
    senimg_name_temp = senimg_name(1:strfind(senimg_name, '.sen_img')); 
    senimg_name_temp = [senimg_name_temp, num2str(start), '_', num2str(finish), '.sen_img'];
    ss = ['now reading in ', senimg_name_temp];
    disp(ss);

    fid = fopen(senimg_name_temp, 'rb');
    sen_temp = fread(fid, inf, 'float');
    sen_img = sen_img + sen_temp;
end

out_name = senimg_name(1:strfind(senimg_name, '.sen_img'));
out_name = [out_name, '_combine_', num2str(num_beds), 'beds_', num2str(ring_per_bed), 'rings_from', num2str(bedStartRing), '.sen_img'];

ss=['output file name ' , out_name]
disp(ss);

fid_out = fopen(out_name, 'w');
fwrite(fid_out, 'float');

sen_img = reshape(sen_img, 239, 239, 679);
slice = reshape(sen_img(:,round(239/2),:), [239, 679]);
imshow(slice, []);
colorbar;

end