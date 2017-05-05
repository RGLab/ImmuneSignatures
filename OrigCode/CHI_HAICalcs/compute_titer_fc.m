function r = compute_titer_fc(initial, final, has_log)

   if ( has_log == 0 ) 
       r.fc = final ./ initial;
   elseif ( has_log ) 
       r.fc = final-initial;
   end
end