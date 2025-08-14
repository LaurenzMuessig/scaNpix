function posOut = quickPosInterpolate(posIn)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%%
arguments
    posIn (:,2) {mustBeNumeric}
end

%%
posOut = posIn;

if ~any(isnan(posIn)); return; end
%
missing_pos   = find(isnan(posIn(:,1)));
ok_pos        = find(~isnan(posIn(:,1)));
for j = 1:2
    posOut(missing_pos, j)                            = interp1(ok_pos, posIn(ok_pos, j), missing_pos, 'linear');
    posOut(missing_pos(missing_pos > max(ok_pos)), j) = posIn( max(ok_pos), j);
    posOut(missing_pos(missing_pos < min(ok_pos)), j) = posIn( min(ok_pos), j);
end



end