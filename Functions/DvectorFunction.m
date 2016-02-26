function [ Dvector ] = DvectorFunction( ref, abs, select_no_atom)
%DvectorFunction Summary of this function goes here
%   Detailed explanation goes here

Dvector = zeros([1 size(ref,3)]);

for j=1:size(ref,3)
Dvector(j) = sum(sum(select_no_atom.*abs.*ref(:,:,j)));
end

Dvector = Dvector';

end