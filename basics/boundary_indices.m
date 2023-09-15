function [i_left, i_right] = boundary_indices(s, x)
% Ќаходит индексы диапазона в списке x, в котором лежит точка s
% —писок x - упор€доченный, значение s лежит внутре

if issorted(x) ~= 0
    if s <= x(1)
        i_right = 1; i_left = i_right;
%         i_left = NaN;
    elseif s > x(end)
        i_left = length(x); i_right = i_left;
%         i_right = NaN;
    else
        for i = 1:(length(x)-1)
            if x(i) < s & x(i+1) >= s
                i_left = i;
                i_right = i+1;
            end
        end
    end
else
    disp('Vector x is not sorted in ascending order.')
end
       
        
        
        