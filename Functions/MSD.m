function [ mean_square ] = MSD( absorption, reference, select_no_atom )
%MSD calculates the mean square deviation between an absorption and
%reference image outside of the signal region defined my the selection
%matrix.

mean_square = sum(sum(...
    select_no_atom .* (squeeze(absorption) - squeeze(reference)).^2));

end

