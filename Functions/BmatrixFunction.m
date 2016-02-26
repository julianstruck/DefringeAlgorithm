function [ Bmatrix ] = BmatrixFunction( ref, select_no_atom)
%Bmatrix This function computes the Bmatrix for the defringe algorithm

Bmatrix = zeros(size(ref,3));

for j=1:size(ref,3)
for i=1:size(ref,3)
Bmatrix(i,j) = sum(sum(select_no_atom.*ref(:,:,i).*ref(:,:,j)));
end
end

end

