function [fMask, mapNorm] = defineMainField( map, options )
% Define main field in pre-barrier trial: this is a couple of lines, but needs 
% to be consistent across several different calling functions.
%
%   [mask] = bvcTrDefineMainField( map, prms )
%
% map is the bsl trial rate map.
% mask is a logical mask for the main field.

%%
arguments
    map {mustBeNumeric} 
    options.fieldThr (1,1) {mustBeNumeric} = 0.5;
    options.mode (1,:) {mustBeMember(options.mode,{'norm2pk','Z'})} = 'Z';
end

%%
if strcmp( options.mode, 'Z' )
    mapNorm = (map-mean(map(:),'omitnan')) ./ std(map(:),'omitnan');
elseif strcmp( options.mode, 'norm2pk' )
    mapNorm = map ./ max(map(:),[],'omitnan');
end

fMask = mapNorm >= options.fieldThr;

end