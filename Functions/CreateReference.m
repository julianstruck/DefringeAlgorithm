function [ matched_reference ] = CreateReference( Cvector , ref)
%CreateReference Summary of this function goes here
%   Detailed explanation goes here

[size1,size2,size3] = size(ref);
matched_reference = zeros(size1,size2);

for i=1:length(Cvector)
    
matched_reference = matched_reference + squeeze(Cvector(i) .* ref(:,:,i));

end

