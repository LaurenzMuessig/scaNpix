function [modInd, gridPropsStrct] = sortModules(dataObj,method,varargin)
%UNTITLED2 Summary of this function goes here
% you need to add the java bits to javaclasspath.txt for it to work


%%
addpath(genpath('D:\Dropbox\matlab\file_exchange\umapAndEppFileExchange\umap\'));
addpath(genpath('D:\Dropbox\matlab\file_exchange\umapAndEppFileExchange\util\'));
addpath(genpath('D:\Dropbox\matlab\file_exchange\umapAndEppFileExchange\epp\'));

%% PRMS
% pairwise
trialInd         = 1;
overlapThresh    = 0.7;
gridnessThresh   = -1;
useLabelGoodOnly = true;
minNspikes       = 300; 
minNCellsModule  = 3;
plotModulesStats = true;
hAx              = "none";
% UMAP
binSzCoarseMap   = 2.5;
radiusExclude    = [15 100];
%
p = inputParser;
addParameter(p,'trialn',trialInd,@isscalar);
addParameter(p,'minspikes',minNspikes,@isscalar);
addParameter(p,'minmod',minNCellsModule,@isscalar);
addParameter(p,'usegood',useLabelGoodOnly,@islogical);
addParameter(p,'overlap',overlapThresh,@isscalar);
addParameter(p,'mingridness',gridnessThresh,@isscalar);
addParameter(p,'binsz',binSzCoarseMap,@isscalar);
addParameter(p,'radius',radiusExclude);
addParameter(p,'plot',plotModulesStats,@islogical);
addParameter(p,'ax',hAx,(@(x) ishghandle(x, 'axes') || isstring(x)));

parse(p,varargin{:});


% mapPrms = scanpix.maps.defaultParamsRateMaps;
% rMapPrms = mapPrms.rate;
switch lower(method) 
    
    case 'pw'
        %
        if isempty(dataObj.maps.sACs{p.Results.trialn})
            scanpix.maps.addMaps(dataObj,'sac',p.Results.trialn);
        end

        % filter for min n spikes
        nSpikes = cellfun(@(x) length(x),dataObj.spikeData.spk_Times{p.Results.trialn});
        
        [tmpGridness,props,propsEllipse] = cellfun(@(x) scanpix.analysis.gridprops(x,'getellgridness',true),dataObj.maps.sACs{p.Results.trialn},'uni',0);
        tmp = cell2mat(tmpGridness);
        gridness = tmp(:,1);
        gridPropsStrct = cell2mat(cellfun(@(x) [x.gridness' x.wavelength' x.orientation' x.offset'],props,'uni',0));
        ellipseFits    = cell2mat(cellfun(@(x) [x.ellOrient x.ellAbScale],propsEllipse,'uni',0));
        %
        cellind = nSpikes > p.Results.minspikes & gridness > p.Results.mingridness; 
        if p.Results.usegood
            cellind = cellind & deblank(dataObj.cell_Label) == 'good';
        end
        % sanity check
        if sum(cellind) < p.Results.minmod
            warning('Less than 4 grid cells found inside dataset');
            modInd = zeros(length(cellind),1);
            modInd(cellind) = 1;
            return
        end
                
        % some init vals for the ellipses
        [cols, rows] = meshgrid(1:length(dataObj.maps.sACs{p.Results.trialn}), 1:length(dataObj.maps.sACs{p.Results.trialn}));
        center = ceil(length(dataObj.maps.sACs{p.Results.trialn})/2);
        
        % Intersection/Union ratio - according to Tocker et al. (2015)
        IUratio = nan(length(dataObj.maps.sACs{p.Results.trialn}),length(dataObj.maps.sACs{p.Results.trialn}));
        for i = 1:length(dataObj.maps.sACs{p.Results.trialn})
            if ~cellind(i) || any(isnan(ellipseFits(i,2:3))); continue; end
            % Create a logical mask of an ellipse with radii 'abScale(1)', 'abScale(2)' and tilt 'orient'
            radiusX = ellipseFits(i,2)/2;
            radiusY = ellipseFits(i,3)/2;
            ellipsePixA = (sin(ellipseFits(i,1)).*(cols - center) + cos(ellipseFits(i,1)).*(rows - center)).^2 ./ radiusX^2 + (cos(ellipseFits(i,1)).*(cols - center) - sin(ellipseFits(i,1)).*(rows - center)) .^2 ./ radiusY^2 <= 1;
            for j = i+1:length(dataObj.maps.sACs{1})
                if ~cellind(j) || any(isnan(ellipseFits(j,2:3))); continue; end
                
                radiusX = ellipseFits(j,2)/2;
                radiusY = ellipseFits(j,3)/2;
                ellipsePixB = (sin(ellipseFits(j,1)).*(cols - center) + cos(ellipseFits(j,1)).*(rows - center)).^2 ./ radiusX^2 + (cos(ellipseFits(j,1)).*(cols - center) - sin(ellipseFits(j,1)).*(rows - center)) .^2 ./ radiusY^2 <= 1; 
                
                IUratio(i,j) = sum(ellipsePixA(:) & ellipsePixB(:)) / sum(ellipsePixA(:) | ellipsePixB(:)); % Intersection / Union 
            end 
        end
        
        % construct graph with the ellipse ratios as edge weights
        adjMat = IUratio; 
        adjMat(isnan(adjMat) | IUratio < p.Results.overlap) = 0;
        G = graph(adjMat,'upper');
        G.Nodes.Name = cellstr(num2str(dataObj.cell_ID(:,1)));
        sG = subgraph(G, degree(G) > 1); % remove all unconnected cells
        
        % split sub graphs if necessary
        if height(sG.Nodes) > 0
            components = conncomp(sG,'output','cell');
            sGraphs = cell(1,length(components));
            for i = 1:length(components)
                if length(components{i}) < p.Results.minmod
                    continue;
                end
                % subfuction
                sGraphs{i} = splitGraph(subgraph(sG,components{i}));
            end
            % format into cell array
            tmpG = sGraphs;
            while any(cellfun(@iscell,tmpG))
                tmpG = [tmpG{cellfun(@iscell,tmpG)} tmpG(~cellfun(@iscell,tmpG))];
            end
            sGraphs = tmpG;
            % generate module index (0=non-grid,1=unsorted grids,2-n=modules)
            modInd = zeros(length(dataObj.cell_ID),1);
            modInd(cellind) = 1;
            for i = 1:length(sGraphs)
                if ~isempty( sGraphs{i})
                    modInd( ismember( G.Nodes.Name, sGraphs{i}.Nodes.Name) ) = i+1;
                end
            end
        else
            % no modules found
            modInd          = zeros(length(dataObj.cell_ID),1);
            modInd(cellind) = 1;
        end
        % 
    
    case 'umap'
          %% Only works with large amount of data and needs some more parameter tuning
%         % pre process
%         % make coarse rate map
%         rMapPrms.binSizeSpat = prms.binSzCoarseMap;
%         rMapPrms.smooth = 'boxcar';
%         rMapPrms.smoothKernel = 1;
%         %
% %         ind = cell2mat(cellfun(@length,dataObj.spikeData.spk_Times{4},'uni',0)) > 150;
%         mapsCoarse = scanpix.maps.makeRateMaps(dataObj.spikeData.spk_Times{4}(286), dataObj.posData.XY{4}, dataObj.spikeData.sampleT{4}, dataObj.trialMetaData(4).ppm, dataObj.posData.speed{4}, rMapPrms );
%         % get sAC
%         sACs = cellfun(@(x) scanpix.analysis.spatialCrosscorr(x, x,'removeMinOverlap',false), mapsCoarse,'uni',0);
%         % remove central peak
%         [x,y] = meshgrid(0.5:size(sACs{1},1)-0.5,0.5:size(sACs{1},2)-0.5);
%         distMap = sqrt( (x - size(x,1)/2).^2 + (y - size(x,2)/2).^2 ) * prms.binSzCoarseMap;
%         indRemove = repmat(distMap < prms.radiusExclude(1) | distMap > prms.radiusExclude(2),1,1,length(sACs));
%         sACs = cat(3,sACs{:});
%         sACs(indRemove|isnan(sACs)) = 0; % is this right? setting to NaN doesn't work for UMAP
%         % reshape (rows = spatial bins, columns = cells)
%         sACs = reshape(sACs,[],size(sACs,3),1);
%         % sACs(all(sACs==0,2),:) = [];
%         % z-score
%         % sACs = (sACs - nanmean(sACs,1)) ./ nanstd(sACs,1);
%         
%         %% umap
%         % this needs the Matlab implementation of UMAP from the file exchange (https://uk.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)
%         [reduction, umap, clusterIdentifiers, extras] = run_umap(sACs','min_dist',0.05,'n_neighbors',5,'metric','cityblock','cluster_detail','very low');
%         % y = tsne(sACs');
%         
%         %% cluster
%         idx = dbscan(reduction,0.8,30);
%         
%         %% ID the grid cell cluster
           
end

if p.Results.plot
    plotModuleStats(dataObj,p.Results.trialn,modInd,gridPropsStrct,IUratio,sGraphs);
end

end

% -------------------------------------------------------------------------------------------------
% --- INLINE FUNCTIONS ----------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
function [sGraphs] = splitGraph(G)
% We split the graph using the Fiedler vector


minNodes = 5;

L = laplacian(G);     % Laplacian matrix of Graph
[V,D] = eig(full(L)); % eigendecomposition of Laplacian
w = V(:,2);           % Fiedler vect

% this indicates that we shouldn't split graph further, a bit arbitrary to
% have a min n of nodes for a splt to occur but the Fiedler method doesn't
% produce great results on a small graph.
% D(2,2)
if D(2,2) > 1 || height(G.Nodes) <= minNodes %
    sGraphs = {G};
    return
end

if find(diff(diag(D)) > 0.05*height(G.Nodes),1,'first') == 2 % somewhat experimental
    if ~any(w<0)
        sGraphs = {subgraph(G,w == 0),subgraph(G,w > 0)};
    elseif ~any(w>0)
        sGraphs = {subgraph(G,w == 0),subgraph(G,w < 0)};
    else
        sGraphs = {subgraph(G,w >= 0),subgraph(G,w < 0)};
    end
else
    sGraphs = {splitGraph(subgraph(G,w >= 0)),splitGraph( subgraph(G,w < 0))}; % BOTH subGraphs need recursive split?
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotModuleStats(dataObj,trialInd,modInd,gridProps,IUratio,graphData)
% plot the results of the module sorting to check that it worked okay 

cols = 'krgbymc';
%
axArr = scanpix.plot.multPlot([5 max([3,length(unique(modInd))-1])],'plotsize',[150 150],'plotsep',[75 40]);
hold(axArr{1,1},'on');
%
withinMod = []; acrossMod = [];
for i = 1:length(unique(modInd))-1
    
    scanpix.plot.mapsMultPlot({dataObj.maps.rate{1}(modInd==i),dataObj.maps.rate{2}(modInd==i)},{'rate'},'cellIDStr',cellstr(num2str(dataObj.cell_ID(modInd==i,1))));
 
    if i > 1
        
        scatter(axArr{1,1},gridProps(modInd==i,3).*180/pi,gridProps(modInd==i,2),'Filled',[cols(i) 'o']);
        
        plot(axArr{1,i},graphData{i-1});
        
        tmp = IUratio(modInd==i,modInd==i);
        withinMod = [withinMod;tmp(~isnan(tmp))];
        %
        tmp = IUratio(modInd==i, modInd~=i & modInd~=0);
        acrossMod = [acrossMod;tmp(~isnan(tmp))];
        %
        histogram(axArr{2,i},gridProps(modInd==i,2),linspace(5,30,13));
        xlabel(axArr{2,i},'grid scale');
        histogram(axArr{3,i},gridProps(modInd==i,3),linspace(0,30*pi/180,16));
        xlabel(axArr{3,i},'grid orientation');
        histogram(axArr{4,i},gridProps(modInd==i,4), linspace(0,30*pi/180,16));
        xlabel(axArr{4,i},'offset');
        %
        imagesc(axArr{5,i},nanmean(cat(3,dataObj.maps.sACs{trialInd}{modInd == i}),3));
        colormap(axArr{5,i},jet);
        axis(axArr{5,i},'square');
    else
        imagesc(axArr{5,1},nanmean(cat(3,dataObj.maps.sACs{trialInd}{modInd == 0}),3)); 
        colormap(axArr{5,1},jet);
        axis(axArr{5,1},'square');
    end
   
end
set(axArr{1,1},'xlim',[-30 30],'ylim',[5 30]);
hold(axArr{1,1},'off');
%
hold(axArr{1,2},'off');
% distribution of across v within module IUratios
histogram(axArr{2,1}, withinMod, linspace(0,1,20),'FaceColor','k');
hold(axArr{2,1},'on');
histogram(axArr{2,1}, acrossMod, linspace(0,1,20),'FaceColor','r');
hold(axArr{2,1},'off');
xlabel(axArr{2,1},'I/U ratio');

end
