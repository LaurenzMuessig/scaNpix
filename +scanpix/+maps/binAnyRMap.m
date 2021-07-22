function [rMapBinned, cMapBinned] = binAnyRMap(rMap,cMap,nSteps,binVals,RGB4nans)
% bin any rate map for plotting and create a corresponding colormap using
% Matlab's in-built colormaps. In particular we want to control the range
% across which we bin, as well as including a color for nan's which are
% normally ignored by colormap().
% Can then plot as imagesc(rMapBinned); colormap(cMapBinned);
% package: scanpix.maps
%
%
% Usage:    
%           [rMapBinned, cMapBinned] = scanpix.maps.binAnyRMap(rMap,cMap,nSteps)
%           [rMapBinned, cMapBinned] = scanpix.maps.binAnyRMap(rMap,cMap,nSteps,RGB4nans)
%           [rMapBinned, cMapBinned] = scanpix.maps.binAnyRMap(rMap,cMap,nSteps,RGB4nans,binVals)
%
%
% Inputs:   rMap     - rate map
%           cMap     - 'string' - name of colormap in matlab
%           nSteps   - how many steps should the colormap for the rate map have?
%           RGB4nans - [R G B]; [1 1 1] as default (white)
%           binVals  - range across which to bin (default: [0 max(rMap(:)])      
%
% Outputs:  rMapBinned - binned rate map
%           cMapBinned - colormap for rate map
%
% LM 2020

if isempty(rMap)
    [rMapBinned, cMapBinned] = deal([]);
    return
end

if nargin < 4
    RGB4nans = [1 1 1];
    binVals  = [0 nanmax(rMap(:))];
elseif nargin == 4
    RGB4nans = [1 1 1];
end

% make colormap
switch lower(cMap)
    case 'hcg' 
        temp = scanpix.maps.highContGrayColMap;
        ind  = round(linspace(1,length(temp),nSteps));
        cMap = colormap(temp(ind,:));   
    case 'poulter'
        cMap = scanpix.maps.cm_Poulter;
    otherwise
        try
            cMap = feval( str2func(cMap), nSteps );
        catch
            error(['''' cMap ''' not yet supported as colormap. Why don''t you add it yourself?']);
        end
end

% bin rMap
rMapBinned              = discretize(rMap,linspace(binVals(1),binVals(2),nSteps+1));
rMapBinned(isnan(rMap)) = 0;

% color map
cMapBinned              = nan(size(cMap,1)+1,size(cMap,2));
cMapBinned(1,:)         = RGB4nans;
cMapBinned(2:end,:)     = cMap;



end

