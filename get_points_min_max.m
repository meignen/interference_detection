function [maxi, mini] = get_points_min_max(sp)
% get_points: extract local minima, local maxima in 2D.

[n,m] = size(sp);

% Mask for Max
maxi =zeros(size(sp));
mini = zeros(size(sp));
for k=2:n-1
 for p=2:m-1
    A = (sp(k,p) - sp(k-1:k+1,p-1:p+1));
    B = A(A >= 0);
    mini(k,p) = (length(B) == 1);
    C = A(A <= 0);
    maxi(k,p) = (length(C) == 1);
 end    
end
end

