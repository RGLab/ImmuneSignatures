function r = robust_med(x)
    r = median( x( ~( isnan(x)|isinf(x) ) ) );
end