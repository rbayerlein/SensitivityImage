close all;

fname = '/home/rbayerlein/Documents/TEMP/CTAC_201.sen_img';
fid = fopen(fname, 'rb');
sen_img = fread(fid, inf, 'float');
sen_img = reshape(sen_img, 239,239,679);

slice = reshape(sen_img(:,round(239/2),:), [239, 679]);


imshow(slice, []);