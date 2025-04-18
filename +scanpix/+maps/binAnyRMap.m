function [rMapBinned, cMapBinned] = binAnyRMap(rMap,options)
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
%           varargin - optional inputs; Name-Value pair(s) for:
%                      'colmap'   - Matlab color map
%                      'nsteps'   - how many steps in colour map), 
%                      'rgb4nans' - RGB triplet for colour for NaNs in map 
%                      'cmapEdge' - [min max] for colormap; default = [0 max(rMap(:)]       
%
% Outputs:  rMapBinned - binned rate map
%           cMapBinned - colormap for rate map
%
% LM 2020

%%
arguments
    rMap {mustBeNumeric}
    options.colmap = 'jet' ;
    options.nsteps (1,1) {mustBeNumeric} = 11;
    options.rgb4nans (1,3) {mustBeNumeric} = [1 1 1];
    options.cmapEdge (1,2) {mustBeNumeric} = [0 max(rMap(:),[],'omitnan')];
end

%%
if isempty(rMap)
    [rMapBinned, cMapBinned] = deal([]);
    return
end

%% make colormap
switch lower(options.colmap)
    case 'hcg' 
        temp   = scanpix.maps.highContGrayColMap;
        ind    = round(linspace(1,length(temp),options.nsteps));
        cMap   = temp(ind,:);   
        nSteps = options.nsteps;
    case 'poulter'
        cMap   = scanpix.maps.cm_Poulter;
        nSteps = size(cMap,1);
    otherwise
        try
            cMap   = feval( str2func(options.colmap), options.nsteps );
            nSteps = options.nsteps;
        catch
            error(['''' options.colmap ''' not yet supported as colormap. Why don''t you add it yourself?']);
        end
end

%% bin rMap
[rMapBinned,edges]          = discretize(rMap,linspace(options.cmapEdge(1),options.cmapEdge(2),nSteps+1));
rMapBinned(rMap<edges(1))   = 1;
rMapBinned(rMap>edges(end)) = max(rMapBinned(:),[],'omitnan');
%
rMapBinned(isnan(rMap))     = 0;

%% color map
cMapBinned              = nan(size(cMap,1)+1,size(cMap,2));
cMapBinned(1,:)         = options.rgb4nans;
cMapBinned(2:end,:)     = cMap;

end

