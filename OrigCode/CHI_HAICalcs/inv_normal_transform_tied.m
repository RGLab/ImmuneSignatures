% 
% for each column of the data, rank them and then use the inverse normal
% transform to enforce the data to be normally distributed in each column
%
function data = inv_normal_transform_tied( in_data )
    data = zeros( size(in_data) );
    rank = zeros( size(in_data,1), 1);
    for i=1:size(in_data,2)
        rank = tiedrank( in_data(:,i) );
        p = rank / ( length(rank) + 1 );
        data(:,i) = norminv( p, 0, 1 );
    end
end