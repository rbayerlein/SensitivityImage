close all;

fname = '/media/rbayerlein/data/recon_data/20200124/sen_img/CTAC_201.combine_1beds_84rings_from85_0%overlap.sen_img';
%fname = '/media/rbayerlein/data/recon_data/20200124/multi_bed/1bed_84rings_from85_120s_4it/lmrecon_explorer_OSEM_f0.intermediate.4';
%fname = '/media/rbayerlein/data/recon_data/20200124/single_bed/0cm_24cm_1200s_4it/lmrecon_explorer_OSEM_f0.intermediate.4_MOD';
fid = fopen(fname, 'rb');
sen_img = fread(fid, inf, 'float');
sen_img = reshape(sen_img, 239,239,679);

slice = reshape(sen_img(:,round(239/2),:), [239, 679]);

imshow(slice, []);

end_of_bed = 1;
for i = 86 : 679
    sum  = 0;
    for x = 1 : 239
        for y = 1 : 239
            sum = sum + sen_img(x,y,i);
        end
    end
    i
    sum
    if sum > 1.0e+12
        sum
        i
        return
    end
end


