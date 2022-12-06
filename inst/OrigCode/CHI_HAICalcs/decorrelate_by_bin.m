function r = decorrelate_by_bin(x, y, top, bot, bins)

   fc = y;
   %diff = final-initial;
   fc_residual = zeros(size(fc))*NaN;
   %diff_residual = zeros(size(fc));
   
   for i=1:length(bins)
       range = bins{i};
       % the range must be non-overlapping
       index = x >= range(1) & x <=range(2);
       %temp = mad(fc(index),1);
    
       fc_residual(index)= (fc(index)-median(fc(index)))/mad(fc(index),1);
       %fc_residual(index)= (fc(index)-mean(fc(index)))/std(fc(index),1);
       %diff_residual(index)= (diff(index)-median(diff(index)))/mad(diff(index),1);
   end
   
   % NaN = 0/0 when the mad = 0
   % inf = k/0 when mad = 0 and k>0
   % -inf = k/0 when mad = 0 and k<0
   % in essence, since there is little variation, we can consider all
   % values to have 0 deviation from the median
   % mask = isnan(fc_residual) | isinf(fc_residual);
   % fc_residual(mask) = 0;
   % mask = isnan(diff_residual) | isinf(diff_residual);
   %diff_residual(mask) = 0;
   
   r.fcd = discretize(fc, top, bot);
   %r.diffd = discretize(diff, top, bot);
   r.fc_res_d = discretize(fc_residual, top, bot);
   %r.diff_res_d = discretize(diff_residual, top, bot);
   
   r.fc = fc;
   %r.diff = diff;
   r.fc_res = fc_residual;
   %r.diff_res = diff_residual;
   
end