function [ Bmatrix ] = Bmatrix( ref, select_no_atom)
%Bmatrix Summary of this function goes here
%   Detailed explanation goes here

B = zeros(size(ref,3));

for j=1:size(ref,3)
for i=1:size(ref,3)
B(i,j) = sum(sum(select_no_atom.*ref(:,:,i).*ref(:,:,j)));
end
end


end

