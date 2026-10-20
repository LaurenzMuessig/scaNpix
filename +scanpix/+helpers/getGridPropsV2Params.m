function gpArgs = getGridPropsV2Params(mapParams)
% getGridPropsV2Params - name-value options for scanpix.analysis.gridprops_v2 from a mapParams struct
% package: scanpix.helpers
%
%       gpArgs = scanpix.helpers.getGridPropsV2Params(obj.mapParams);
%       [g, Props] = scanpix.analysis.gridprops_v2(sac, gpArgs{:});
%
% Uses mapParams.gridPropsV2; objects saved before this field existed fall back to the
% defaults in scanpix.maps.defaultParamsRateMaps.
%
% LM / Claude 2026

if nargin > 0 && isstruct(mapParams) && isfield(mapParams, 'gridPropsV2') && ~isempty(mapParams.gridPropsV2)
    prms = mapParams.gridPropsV2;
else
    def  = scanpix.maps.defaultParamsRateMaps;
    prms = def.gridPropsV2;
end
gpArgs = namedargs2cell(prms);
end
