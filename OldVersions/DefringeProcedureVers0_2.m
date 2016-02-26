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


%%% Preparatory part of the program

%% define file paths
addpath('/Users/Julian/Documents/MIT/MatlabPrograms/MatlabFunctionsJulian')
addpath(fullfile(pwd,'Functions'))
raw_data_path = fullfile(pwd,'RawData');
processed_data_path = fullfile(pwd,'ProcessedData');


%% load raw data

[RawImgName,RawImgPathName,RawData,LoadSuccess] = LoadFitsSeries(raw_data_path);

%% create OD images (no saturation effects included)
OD_images_unmodified = squeeze(Raw2ODseries( RawData(:,:,:,:) ));


%% define picture frame (noatoms and atoms)

% If necessary crop the whole image
OD_image_crop_select  = squeeze(sum(OD_images_unmodified,3));
[crop_select, rect_crop] = imcrop(OD_image_crop_select);
rect_crop=round(rect_crop);
crop_region_y = rect_crop(1):(rect_crop(1) + rect_crop(3) - 1);
crop_region_x = rect_crop(2):(rect_crop(2) + rect_crop(4) - 1);

crop_raw_data = RawData(crop_region_x,crop_region_y,:,:);

% select region with absorption signal
[atom_select, rect_atom] = imcrop(crop_select);
rect_atom=round(rect_atom);
atom_crop_y = rect_atom(1):(rect_atom(1) + rect_atom(3) - 1);
atom_crop_x = rect_atom(2):(rect_atom(2) + rect_atom(4) - 1);

%% define no atom selection matrix
select_no_atom = ones(length(crop_region_x),length(crop_region_y));
select_no_atom(atom_crop_x,atom_crop_y)=0;

%% correct the raw absorption and reference images with their dark frames
% Absorption image
absorption = squeeze(crop_raw_data(:,:,1,:));
    
% No atoms reference
reference = squeeze(crop_raw_data(:,:,2,:));
   
% Dark frame for absoprtion image
dark_abs = squeeze(crop_raw_data(:,:,3,:));
     
% Dark frame for reference image
dark_ref = squeeze(crop_raw_data(:,:,4,:));

% subtract dark image
cor_abs = absorption - dark_abs;
cor_abs(cor_abs<0.1)=0.1;
     
cor_ref = reference - dark_ref;
cor_ref(cor_ref<0.1)=0.1;

%% Main part of the program

% create Bmatrix
Bmatrix = BmatrixFunction( cor_ref, select_no_atom);

matched_reference_vec = zeros([size(cor_ref,1),size(cor_ref,2),length(RawImgName)]);
OD_image = matched_reference_vec;
coefficient_vec = zeros(size(cor_ref,3));

for i=1:length(RawImgName)
% create Dvector
Dvector = DvectorFunction( cor_ref, cor_abs(:,:,i), select_no_atom);

% find coefficients for optimal reference image
Cvector = linsolve(Bmatrix,Dvector);

matched_reference = CreateReference( Cvector , cor_ref);

matched_reference_vec(:,:,i) = matched_reference;
coefficient_vec(:,i) = Cvector;
OD_image(:,:,i) = - log( cor_abs(:,:,i) ./ matched_reference_vec(:,:,i) );

end

save('defringed_images.mat','OD_image')

imagesc(OD_image(:,:,1))
imagesc(squeeze(Raw2OD( crop_raw_data(:,:,:,1)) ))
