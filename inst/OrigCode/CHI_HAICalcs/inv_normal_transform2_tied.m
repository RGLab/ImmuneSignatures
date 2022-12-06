% 
% this version ignores NaNs and return a inverse normal for the non-NaN
% values and keeps the NaNs intact
%
function data = inv_normal_transform2_tied( in_data )
    data = zeros( size(in_data) );
    for i=1:size(in_data,2)
        mask = ~isnan(in_data(:,i));
        data(~mask,i) = NaN;
        data(mask,i) = inv_normal_transform_tied( in_data(mask,i) );
    end
end