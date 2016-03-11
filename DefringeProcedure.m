%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defringe program
%
% This program generates an ideal reference image for a given absorption
% image based on a large basis set of available reference images. This is
% done by minimizing the least square between a given absorption image and 
% a linear combination of all available reference images.
%
% The actual procedure is based on PRA 82, 061606(R) (2010).
% Outline of the procedure:
% Optimized reference image Q_i for absorption image A_i and set of
% reference image {R_k}
% Q_i(x,y) = sum_k c_k*R(x,y)_k 
% Finding the right coefficients c_k is the main task of this program.
% The right set of coefficient is determined by method of least square
% sum_(x,y) m(x,y)(A_i(x,y) - Q_i(x,y))^2, 
% where m(x,y) crops out the region of the image with the atoms (m=0).
% Setting the partial derivatives of the above with respect to the 
% coefficients to zero results in
% sum_k c_k B_(j,k) = D_j,
% where B_(j,k) = sum_(x,y) m(x,y)*R_j(x,y)*R_k(x,y) and 
% D_j=sum_(x,y) m(x,y)*R_j(x,y)*A(x,y).
% This linear system of equations can be readliy solved using linsolve to
% obtain the coefficients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

%%% Preparatory part of the program

%% define file paths
addpath('/Users/Julian/Documents/MIT/MatlabPrograms/MatlabFunctionsJulian')
addpath(fullfile(pwd,'Functions'))
raw_data_path = fullfile(pwd,'TestRawData');
processed_data_path = fullfile(pwd,'ProcessedData');


%% load raw data

[RawImgName,RawImgPathName,RawData,LoadSuccess] = LoadFitsSeries(raw_data_path);

%% create OD images (no saturation effects included)
OD_images_unmodified = squeeze(Raw2ODseries( RawData(:,:,:,:) ));


%% define picture frame (noatoms and atoms)
% sum up all images to see where atoms are;
OD_image_crop_select  = squeeze(sum(OD_images_unmodified,3))/size(RawData,4);
figure
imagesc(OD_image_crop_select,[-0.3 0.2])
axis image;

% If necessary crop the whole image (precrop)
% Construct a question dialog
choice = questdlg('Would you like to precrop the image?', ...
	'Precrop Menu', ...
	'Yes','No','No');
% Handle response
switch choice
    case 'Yes'
        close
        [crop_select, rect_crop] = imcrop(OD_image_crop_select);
        rect_crop=round(rect_crop);
        crop_region_y = rect_crop(1):(rect_crop(1) + rect_crop(3) - 1);
        crop_region_x = rect_crop(2):(rect_crop(2) + rect_crop(4) - 1);
        crop_raw_data = RawData(crop_region_x,crop_region_y,:,:);
    case 'No'
        close
        crop_raw_data = RawData(:,:,:,:);
        crop_select = OD_image_crop_select;
        crop_region_y = 1:size(crop_raw_data,1);
        crop_region_x = 1:size(crop_raw_data,2);
end

% select region with absorption signal
[atom_select, rect_atom] = imcrop(crop_select);
rect_atom=round(rect_atom);
atom_crop_y = rect_atom(1):(rect_atom(1) + rect_atom(3) - 1);
atom_crop_x = rect_atom(2):(rect_atom(2) + rect_atom(4) - 1);

%% define no atom selection matrix
select_no_atom = ones(length(crop_region_x),length(crop_region_y));
select_no_atom(atom_crop_x,atom_crop_y)=0;

%% correct the raw absorption and reference images with their dark frames
[OD_images_cropped,cor_abs,cor_ref] = Raw2ODseries(crop_raw_data);

%% Main part of the program

% create Bmatrix
Bmatrix = BmatrixFunction( cor_ref, select_no_atom);

% preallocate for loop
matched_reference_vec = zeros([size(cor_ref,1),size(cor_ref,2),length(RawImgName)]);
OD_defringed = matched_reference_vec;
coefficient_vec = zeros(size(cor_ref,3));
FitsArray = zeros([size(cor_ref,1),size(cor_ref,2),2,length(RawImgName)]);
mean_square_matched = zeros([1 length(RawImgName)]);
mean_square_original = zeros([1 length(RawImgName)]);

h = waitbar(0,'Generating optimal reference images');

for i=1:length(RawImgName)
% create Dvector
Dvector = DvectorFunction( cor_ref, cor_abs(:,:,i), select_no_atom);

% find coefficients for optimal reference image
Cvector = linsolve(Bmatrix,Dvector);

matched_reference = CreateReference( Cvector , cor_ref);

matched_reference_vec(:,:,i) = matched_reference;
coefficient_vec(:,i) = Cvector;
OD_defringed(:,:,i) = - log( cor_abs(:,:,i) ./ matched_reference_vec(:,:,i) );

% quality check: compute mean square deviation for matched reference and
% original reference image
mean_square_matched(i) = MSD( cor_abs(:,:,i), matched_reference_vec(:,:,i), select_no_atom );
mean_square_original(i) = MSD( cor_abs(:,:,i), cor_ref(:,:,i), select_no_atom );

% create image array for .fits file
FitsArray(:,:,1,i) = cor_abs(:,:,i);
FitsArray(:,:,2,i) = matched_reference_vec(:,:,i);

waitbar(i/length(RawImgName))
end
close(h)

%% Save defringed images as .fits
save('defringed_images.mat','OD_defringed')

for i=1:length(RawImgName)
    [pathstr,DateName,ext] = fileparts(RawImgName{i});
    filename = [DateName,'_defringed'];
    fitswrite(FitsArray(:,:,:,i),fullfile(processed_data_path,filename))
end

%% Show one example image
figure(1)
imagesc(OD_defringed(:,:,1),[0 0.5])
%imagesc(squeeze(Raw2OD( crop_raw_data(:,:,:,1)) ))
figure(2)
plot(mean_square_matched./mean_square_original)
