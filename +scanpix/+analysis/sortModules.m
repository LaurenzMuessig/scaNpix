function [modInd, gridPropsStrct] = sortModules(dataObj,method,varargin)
%UNTITLED2 Summary of this function goes here
% you need to add the java bits to javaclasspath.txt for it to work


%%
addpath(genpath('D:\Dropbox\matlab\file_exchange\umapFileExchange\umap\'));
addpath(genpath('D:\Dropbox\matlab\file_exchange\umapFileExchange\util\'));
addpath(genpath('D:\Dropbox\matlab\file_exchange\umapFileExchange\epp\'));

%% PRMS
% pairwise
trialInd            = 1;
overlapThresh       = [0.6 0.7];
propLowerThanThresh = 0.1;
gridnessThresh      = 0.5;
useLabelGoodOnly    = true;
minNspikes          = 300; 
minNCellsModule     = 3;
plotModulesStats    = false;
hAx                 = "none";
% UMAP
binSzCoarseMap      = 5;
radiusExclude       = [15 100];
min_dist            = 0.05;
nNeighbours         = 5;
%
p = inputParser;
addParameter(p,'trialn',trialInd,@isscalar);
addParameter(p,'minspikes',minNspikes,@isscalar);
addParameter(p,'minmod',minNCellsModule,@isscalar);
addParameter(p,'usegood',useLabelGoodOnly,@islogical);
addParameter(p,'overlap',overlapThresh,@isscalar);
addParameter(p,'proplowthr',propLowerThanThresh,@isscalar);
addParameter(p,'mingridness',gridnessThresh,@isscalar);
addParameter(p,'binsz',binSzCoarseMap,@isscalar);
addParameter(p,'radius',radiusExclude);
addParameter(p,'mindist',min_dist,@isscalar);
addParameter(p,'nn',nNeighbours,@isscalar);
addParameter(p,'plot',plotModulesStats,@islogical);
addParameter(p,'ax',hAx,(@(x) ishghandle(x, 'axes') || isstring(x)));

parse(p,varargin{:});

% mapPrms = scanpix.maps.defaultParamsRateMaps;
% rMapPrms = mapPrms.rate;

if isempty(dataObj.maps.sACs{p.Results.trialn})
    scanpix.maps.addMaps(dataObj,'sac',p.Results.trialn,scanpix.maps.defaultParamsRateMaps);
end

%% filter cells and get grid props
% filter for min n spikes
nSpikes = cellfun(@(x) length(x),dataObj.spikeData.spk_Times{p.Results.trialn});

[tmpGridness,props] = cellfun(@(x) scanpix.analysis.gridprops(x,'getellgridness',true),dataObj.maps.sACs{p.Results.trialn},'uni',0);
tmp = cell2mat(tmpGridness);
gridness = max(tmp,[],2,'omitnan');
gridPropsStrct = cell2mat(cellfun(@(x) [x.gridness(1,1)' x.wavelength(1,1)' x.orientation(1,1) x.offset(1,1)'],props,'uni',0));
ellipseFits    = cell2mat(cellfun(@(x) [x.ellOrient(1,2) x.ellAbScale(1,3:4)],props,'uni',0));
%
cellind = nSpikes > p.Results.minspikes & gridness > p.Results.mingridness;
if p.Results.usegood
    cellind = cellind & deblank(dataObj.cell_Label) == 'good';
end
% sanity check
if sum(cellind) < p.Results.minmod
    warning('scaNpix::analysis::sortModules: Less than %i grid cells found inside dataset',p.Results.minmod);
    modInd = zeros(length(cellind),1);
    modInd(cellind) = 1;
    return
end
% get Intersection/Union ratio
IUratio = scanpix.analysis.computeIURatio(dataObj,p.Results.trialn,ellipseFits,cellind);

%% assign modules
switch lower(method) 
    
    case 'pw'
        %
        % construct graph with the ellipse ratios as edge weights
        adjMat = IUratio; 
        adjMat(isnan(adjMat) | IUratio < p.Results.overlap(2)) = 0;
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
                tempG      = subgraph(sG,components{i});
                ind        = ismember(G.Nodes.Name,tempG.Nodes.Name);
                sGraphs{i} = splitGraph(tempG,IUratio(ind,ind),p);
            end
            % format into cell array
            tmpG = sGraphs;
            while any(cellfun(@iscell,tmpG))
                tmpG = [tmpG{cellfun(@iscell,tmpG)} tmpG(~cellfun(@iscell,tmpG))];
            end
            sGraphs         = tmpG;
            % generate module index (0=non-grid,1=unsorted grids,2-n=modules)
            modInd          = zeros(length(dataObj.cell_ID),1);
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
            warning('scaNpix::analysis::sortModules: No modules found in dataset');
            return 
        end
        % 
    
    case 'umap'
        %% doesn't work too well with small datasets; post umap clustering needs work 
        % pre process
        % make coarse rate map
        rMapPrms.binSizeSpat  = p.Results.binsz;
        rMapPrms.smooth       = 'boxcar';
        rMapPrms.smoothKernel = 1;
        rMapPrms.showWaitBar  = true;
        %
        mapsCoarse        = scanpix.maps.makeRateMaps(dataObj.spikeData.spk_Times{1}, dataObj.posData.XY{1}, dataObj.spikeData.sampleT{1}, dataObj.trialMetaData(1).ppm, dataObj.posData.speed{1}, rMapPrms );
        % get sAC
        sACs              = cellfun(@(x) scanpix.analysis.spatialCrosscorr(x, x,'removeMinOverlap',false), mapsCoarse,'uni',0);
        % remove central peak
        [x,y]             = meshgrid(0.5:size(sACs{1},1)-0.5,0.5:size(sACs{1},2)-0.5);
        distMap           = sqrt( (x - size(x,1)/2).^2 + (y - size(x,2)/2).^2 ) * p.Results.binsz;
        indRemove         = repmat(distMap < p.Results.radius(1) | distMap > p.Results.radius(2),1,1,length(sACs));
        sACs              = cat(3,sACs{:});
        sACs(indRemove)   = 0; % having NaNs doesn't work for UMAP
        % reshape (rows = spatial bins, columns = cells)
        sACs              = reshape(sACs,[],size(sACs,3),1);
        % sACs(all(sACs==0,2),:) = [];
        % z-score
        sACs              = (sACs - mean(sACs,1,'omitnan')) ./ std(sACs,1,'omitnan');
        sACs(isnan(sACs)) = 0; 
        
        %% umap
        % this needs the Matlab implementation of UMAP from the file exchange (https://uk.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)
        [reduction, umap, clusterIdentifiers, extras] = run_umap(sACs','min_dist',p.Results.mindist,'n_neighbors',p.Results.nn,'metric','cityblock','cluster_detail','very low','verbose','text');
%         % y = tsne(sACs');
        cluIDs = unique(clusterIdentifiers);
        cluIDs(cluIDs==0) = [];
        %
        meanOverlap = nan(length(cluIDs), 1);
        for i = cluIDs
            meanOverlap(i) = mean(mean(IUratio( clusterIdentifiers' == i & cellind, clusterIdentifiers'== i & cellind ),'omitnan'),'omitnan');
        end
        [~, sortInd] = sort(meanOverlap);
        %
        modInd = zeros(length(dataObj.cell_ID), 1);
        c = 1;
        for i = sortInd'
            modInd( clusterIdentifiers' == i & cellind ) = c;
            c = c+1;
        end
        modInd(clusterIdentifiers == 0) = 0;
        %
        sGraphs = [];
%         %% cluster
%         idx = dbscan(reduction,0.8,30);
%         
%         %% ID the grid cell cluster
end

if p.Results.plot
    plotModuleStats(dataObj, p.Results.trialn, modInd, gridPropsStrct, IUratio, sGraphs, p.Results.overlap(2));
end

end

% -------------------------------------------------------------------------------------------------
% --- INLINE FUNCTIONS ----------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
function [sGraphs] = splitGraph(G, IURatio, params)
% We split the graph using the Fiedler method

if height(G.Nodes) <= params.Results.minmod %
    sGraphs = {G};
    return
end

L     = laplacian(G); % Laplacian matrix of Graph
[V,~] = eig(full(L)); % eigendecomposition of Laplacian
w     = V(:,2);       % Fiedler vect

% decide if we need to split the graph further (recursively) - this is
% based on proportion of bad overlaps in distribution 
if sum(IURatio(:) < params.Results.overlap(1)) / sum(~isnan(IURatio(:))) > params.Results.proplowthr
    mod1_ind = ismember(G.Nodes.Name,G.Nodes.Name(w >= 0));
    mod2_ind = ismember(G.Nodes.Name,G.Nodes.Name(w < 0));
    sGraphs  = {splitGraph(subgraph(G,w >= 0),IURatio(mod1_ind,mod1_ind),params),splitGraph(subgraph(G,w < 0),IURatio(mod2_ind,mod2_ind),params)};
else
%     sGraphs = {subgraph(G,degree(G) > ceil(prctile(degree(G),1)))};
    sGraphs  = {G};
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotModuleStats(dataObj, trialInd, modInd, gridProps, IUratio, graphData, thr)
% plot the results of the module sorting to check that it did a sensible job 

cols = 'krgbymc';

%
withinMod = []; acrossMod = [];
modInData = unique(modInd);
modInData = modInData(modInData~=0)';
if ~any(modInData==1)
    modInData = [1, modInData];
end
%
axArr = scanpix.plot.multPlot([5 length(modInData)],'plotsize',[150 150],'plotsep',[75 40]);
hold(axArr{1,1},'on');
%
pltInd = 1;
for i = modInData
    %%%%% THIS IS WRONG AS IT PLOTS THE TRIAL TWICE NEEDS CHANGING %%%
    scanpix.plot.mapsMultPlot({dataObj.maps.rate{trialInd}(modInd==i),dataObj.maps.rate{trialInd}(modInd==i)},{'rate'},'cellIDStr',cellstr(num2str(dataObj.cell_ID(modInd==i,1))));
 
    if i > 1
        
        tmpOrient = gridProps(modInd==i,3).*180/pi;
        tmpOrient( gridProps(modInd==i,3).*180/pi<0) = tmpOrient( gridProps(modInd==i,3).*180/pi<0) + 30;
        tmpOrient( gridProps(modInd==i,3).*180/pi>0) = tmpOrient( gridProps(modInd==i,3).*180/pi>0) - 30;
        scatter(axArr{1,1},tmpOrient,gridProps(modInd==i,2)*2.5,'Filled',[cols(i) 'o']);
        %
        if ~isempty(graphData)
            plot(axArr{1,pltInd},graphData{pltInd-1});
        end
        
        tmp = IUratio(modInd==i,modInd==i);
        withinMod = [withinMod;tmp(~isnan(tmp))];
        %
        tmp = IUratio(modInd==i, modInd~=i & modInd>1);
        acrossMod = [acrossMod;tmp(~isnan(tmp))];
        %
        histogram(axArr{2,pltInd},gridProps(modInd==i,2),linspace(5,30,13));
        xlabel(axArr{2,pltInd},'grid scale');
        histogram(axArr{3,pltInd},gridProps(modInd==i,3),linspace(-180*pi/180,180*pi/180,32));
        xlabel(axArr{3,pltInd},'grid orientation');
        histogram(axArr{4,pltInd},gridProps(modInd==i,4), linspace(0,30*pi/180,16));
        xlabel(axArr{4,pltInd},'offset');
        %
        imagesc(axArr{5,pltInd},nanmean(cat(3,dataObj.maps.sACs{trialInd}{modInd == i}),3));
        colormap(axArr{5,pltInd},jet);
        axis(axArr{5,pltInd},'square');
    else
        tmpOrient = gridProps(modInd==i,3).*180/pi;
        tmpOrient( gridProps(modInd==i,3).*180/pi<0) = tmpOrient( gridProps(modInd==i,3).*180/pi<0) + 30;
        tmpOrient( gridProps(modInd==i,3).*180/pi>0) = tmpOrient( gridProps(modInd==i,3).*180/pi>0) - 30;
        scatter(axArr{1,1},tmpOrient,gridProps(modInd==i,2)*2.5,'Filled','ko');
        %
        imagesc(axArr{5,1},mean(cat(3,dataObj.maps.sACs{trialInd}{modInd == 0}),3),'omitnan'); 
        colormap(axArr{5,1},jet);
        axis(axArr{5,1},'square');
    end
    pltInd = pltInd + 1;
end
plot(axArr{1,1},[-30 -30; 30 30]',[0 62.5; 0 62.5]','k--');
set(axArr{1,1},'xlim',[-180 180],'ylim',[12.5 62.5],'xtick',-180:10:180,'xticklabel',{'0','10','20','30/-30','-20','-10','0'},'XTickLabelRotation',45);
hold(axArr{1,1},'off');
ylabel(axArr{1,1},'grid scale (cm)'); xlabel(axArr{1,1},'grid orientation (deg)');
%
hold(axArr{1,2},'off');
% distribution of across v within module IUratios
histogram(axArr{2,1}, withinMod, linspace(0,1,20),'FaceColor','k');
hold(axArr{2,1},'on');
histogram(axArr{2,1}, acrossMod, linspace(0,1,20),'FaceColor','r');
plot(axArr{2,1},[thr thr]',[0 max(ylim(axArr{2,1}))]','k--');
hold(axArr{2,1},'off');
xlabel(axArr{2,1},'I/U ratio');

end


% if ~any(w<0)
%     sGraphs = {subgraph(G,w == 0),subgraph(G,w > 0)};
% elseif ~any(w>0)
%     sGraphs = {subgraph(G,w == 0),subgraph(G,w < 0)};
% else
%     sGraphs = {subgraph(G,all(w >= 0,2)),subgraph(G,all(w < 0,2))};
% % end
% ind = ismember(G.Nodes.Name,G.Nodes.Name(all(w >= 0,2)));
% splitGraph(sGraphs{1},IURatio(ind,ind));


% if D(2,2) > 1 || height(G.Nodes) <= minNodes %
% % if height(G.Nodes) <= minNodes %
%     sGraphs = {G};
%     return
% end
% 
% if find(diff(diag(D)) > 0.05*height(G.Nodes),1,'first') == 2 % somewhat experimental
%     if ~any(w<0)
%         sGraphs = {subgraph(G,w == 0),subgraph(G,w > 0)};
%     elseif ~any(w>0)
%         sGraphs = {subgraph(G,w == 0),subgraph(G,w < 0)};
%     else
%         sGraphs = {subgraph(G,w >= 0),subgraph(G,w < 0)};
%     end
% else
%     sGraphs = {splitGraph(subgraph(G,w >= 0)),splitGraph( subgraph(G,w < 0))}; % BOTH subGraphs need recursive split?
% end
