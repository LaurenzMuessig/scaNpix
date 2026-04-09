function mapNorm = make_zMap( map, mode )
% Define main field in pre-barrier trial: this is a couple of lines, but needs 
% to be consistent across several different calling functions.
%
%

%%
arguments
    map {mustBeNumeric} 
    mode (1,:) {mustBeMember(mode,{'norm2pk','Z'})} = 'Z';
end

%%
if isempty(map)
    mapNorm = [];
    return
end
%
if strcmp( mode, 'Z' )
    mapNorm = (map-mean(map(:),'omitnan')) ./ std(map(:),[],'omitnan');
elseif strcmp( mode, 'norm2pk' )
    mapNorm = map ./ max(map(:),[],'omitnan');
end

end