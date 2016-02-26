function [ OD_image , correct_absorption , correct_reference] =...
    Raw2ODseries( RawData )
% Raw2OD This function converts a series of raw files into a series of OD 
% images. Input convention: RawData(x,y,ImageType,SeriesNumber).
   
% check size of ImageType and build case structure for 3 and 4
number_of_images = size(RawData,3);

switch number_of_images
    case 3
     % Absorption image
     absorption = squeeze(RawData(:,:,1,:));
    
     % No atoms reference
     reference = squeeze(RawData(:,:,2,:));
   
     % Dark frame
     dark = squeeze(RawData(:,:,3,:));

     % subtract dark image
     correct_absorption = absorption - dark;
     correct_absorption(correct_absorption<0.1)=0.1;
     
     correct_reference = reference - dark;
     correct_reference(correct_reference<0.1)=0.1;
     
    case 4
     % Absorption image
     absorption = squeeze(RawData(:,:,1,:));
    
     % No atoms reference
     reference = squeeze(RawData(:,:,2,:));
   
     % Dark frame for absoprtion image
     dark_abs = squeeze(RawData(:,:,3,:));
     
     % Dark frame for reference image
     dark_ref = squeeze(RawData(:,:,4,:));

     % subtract dark image
     correct_absorption = absorption - dark_abs;
     correct_absorption(correct_absorption<0.1)=0.1;
     
     correct_reference = reference - dark_ref;
     correct_reference(correct_reference<0.1)=0.1;
        
    otherwise
     disp('number of images per shot not supported')
        
end
     % calculate Atom number (not taking into saturation effects)
     OD_image = - log( correct_absorption ./ correct_reference );

end

