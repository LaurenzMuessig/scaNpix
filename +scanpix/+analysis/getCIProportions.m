function ci95 = getCIProportions(n,nTotal)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% from Zarr p.??
q    = n ./ nTotal;
qHat = 1 - q;
ci95 = 2 .* (1.96 .* sqrt( ( q.*qHat ) ./ nTotal )); % two-tailed

end