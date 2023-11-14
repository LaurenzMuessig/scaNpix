function [rMapBinned, cMapBinned] = binAnyRMap(rMap,varargin)
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
defaultColMap    = 'jet';
defaultNSteps    = 11;
defaultRGB4nans  = [1 1 1]; % unvisited bins = white
defaultCMapEdge  = [];

p = inputParser;
addParameter(p,'colmap',defaultColMap,@ischar);
addParameter(p,'nsteps',defaultNSteps,@isscalar); 
addParameter(p,'rgb4nans',defaultRGB4nans);
addParameter(p,'cmapEdge',defaultCMapEdge);
parse(p,varargin{:});

%%
if isempty(rMap)
    [rMapBinned, cMapBinned] = deal([]);
    return
end

% make colormap
switch lower(p.Results.colmap)
    case 'hcg' 
        temp   = scanpix.maps.highContGrayColMap;
        ind    = round(linspace(1,length(temp),p.Results.nsteps));
        cMap   = temp(ind,:);   
        nSteps = p.Results.nsteps;
    case 'poulter'
        cMap   = scanpix.maps.cm_Poulter;
        nSteps = size(cMap,1);
    otherwise
        try
            cMap   = feval( str2func(p.Results.colmap), p.Results.nsteps );
            nSteps = p.Results.nsteps;
        catch
            error(['''' p.Results.colmap ''' not yet supported as colormap. Why don''t you add it yourself?']);
        end
end

% bin rMap
if isempty(p.Results.cmapEdge)
    rMapBinned                  = discretize(rMap,linspace(0,nanmax(rMap(:)),nSteps+1));
else
    [rMapBinned,edges]          = discretize(rMap,linspace(p.Results.cmapEdge(1),p.Results.cmapEdge(2),nSteps+1));
    rMapBinned(rMap<edges(1))   = 1;
    rMapBinned(rMap>edges(end)) = nanmax(rMapBinned(:));
end
%
rMapBinned(isnan(rMap)) = 0;

% color map
cMapBinned              = nan(size(cMap,1)+1,size(cMap,2));
cMapBinned(1,:)         = p.Results.rgb4nans;
cMapBinned(2:end,:)     = cMap;



end

