

clear;




%%% Prepatory part of the program

%% define file paths
addpath('/Users/Julian/Documents/MIT/MatlabPrograms/MatlabFunctionsJulian')
addpath(fullfile(pwd,'Functions'))
raw_data_path = fullfile(pwd,'RawData');
processed_data_path = fullfile(pwd,'ProcessedData');


%% load raw data

[RawImgName,RawImgPathName,RawData,LoadSuccess] = LoadFitsSeries(raw_data_path);

%% define picture frame (noatoms and atoms)
% for now fixed can be replaced with gui input
%OD_image  = Raw2OD( RawData(:,:,:,1) );
%imagesc(OD_image)

crop_region_x = 1:512;
crop_region_y = 10:256; 

crop_raw_data = RawData(crop_region_x,crop_region_y,:,:);

atom_crop_x = 187:354;
atom_crop_y = 53:212;

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
