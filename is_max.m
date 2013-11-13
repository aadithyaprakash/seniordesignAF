function [ result ] = is_max( original, right, left )
% Check to see if the point is greater than its neighbors
    if(original > right && original > left)
        result = 1;
    else
        result = 0;
    end
end

