function r = mad_standardize( x, mad_flag )
    if ~exist('mad_flag', 'var')
        mad_flag = 1; % as in the original function (by John)
    end
    r=(x-robust_med(x))/mad(x,mad_flag); % change flag to 0 (for mean AD) to avoid NaN if more than half samples not changed at all
end
