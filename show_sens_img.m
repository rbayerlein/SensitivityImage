close all;

fname = '/home/rbayerlein/Documents/Projects/20210515_Multi_Bed_Imaging/sens_data_temp/CTAC_201.comb.sen_img';
fid = fopen(fname, 'rb');
sen_img = fread(fid, inf, 'float');
sen_img = reshape(sen_img, 239,239,679);

slice = reshape(sen_img(:,round(239/2),:), [239, 679]);


imshow(slice, []);