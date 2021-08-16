function [mu_map, img_pars] = make_mumap(AC_imgpath)

%ip = genpath('../../image'); 
%addpath(ip); 

[img, CT_AC_dcminfo] = dcm2raw_fast(AC_imgpath, 2); 

%rmpath(ip); 

img_pars = struct; 

img_pars.img_size = size(img);  
img_pars.vox_size = [double(CT_AC_dcminfo.PixelSpacing(1)) double(CT_AC_dcminfo.PixelSpacing(1)) double(CT_AC_dcminfo.SliceThickness)]
img_pars.CT_kVp = double(CT_AC_dcminfo.KVP)

while img_pars.img_size(3) > 828.5
	img = img(:,:,1:2:end) + img(:,:,2:2:end); 
	img = img ./ 2; 
	img_pars.vox_size(3) = img_pars.vox_size(3) * 2;
	img_pars.img_size = size(img)
end

while img_pars.img_size(1) > 256.5

	img = img(1:2:end,:,:) + img(2:2:end,:,:); 
	img = img(:,1:2:end,:) + img(:,2:2:end,:); 
	img = img ./ 4; 
	img_pars.img_size = size(img); 

	img_pars.vox_size(1) = img_pars.vox_size(1)*2; 
	img_pars.vox_size(2) = img_pars.vox_size(2)*2;
end

% shift by -3 mm in z, -2.5 mm in y (fixed offset between UIH CT and UCD PET recon
if abs(img_pars.vox_size(3)  - 3) < 0.5 && abs(img_pars.vox_size(1)  - 2.5) < 0.5
	img(:,:,1:(end-1)) = img(:,:, 2:end); 
	img(:,1:(end-1),:) = img(:,2:end,:); 
end


% rotate the image to match UCD reconstruction
%for zz = 1:img_pars.img_size(3)
%	img(:,:,zz) = rot90(img(:,:,zz), 2); 
%end

img = imgaussfilt3(img, [3/(2.355*img_pars.vox_size(1)) 3/(2.355*img_pars.vox_size(1)) 3/(2.355*img_pars.vox_size(3))]); 
img(img<-800) = -1000; 

% bi-linear 511 keV interpolation
b = 4.71e-3; 
a = 5.1e-6; 

if img_pars.CT_kVp == 80
	b = 6.26e-3; 
	a = 3.64e-6; 
elseif img_pars.CT_kVp  == 120
	b = 4.71e-3; 
	a = 5.1e-6; 
elseif img_pars.CT_kVp == 140
	b = 4.08e-3; 
	a = 5.64e-6; 
else
	errordlg('Invalid CT kVp'); 
	return; 
end

mu_map = img; 
mu_map(img<40) = (img(img < 40) + 1000) .* (9.6e-6); 
mu_map(img>=40) = b + (a .* (img(img >= 40)+1000)); 








