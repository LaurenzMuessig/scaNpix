function ResT = obj2table(obj,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% TODO %
% deal with reordering trials in output T to standard seq

%
scores = {'SI','RV','gridness','borderScore','intraStab'};

%% Params
rowFormat    = 'cell';
maxNTrial    = [];
minNSpikes   = [];
trialIndex   = [];
addMaps      = true;
addScores    = true(1,length(scores));

%
p = inputParser;
addOptional(p, 'rowformat', rowFormat,    ( @(x) mustBeMember(x,{'cell','dataset'}) ));
addParameter(p,'maxn',      maxNTrial,    ( @(x) isscalar(x) || isempty(x) ) );
addParameter(p,'minspikes', minNSpikes,   ( @(x) isscalar(x) || isempty(x) ) );
addParameter(p,'tindex',    trialIndex,   ( @(x) isempty(x) || ischar(x) || isstring(x) ) );
addParameter(p,'addmaps',   addMaps,      @islogical );
addParameter(p,'scores',    addScores,    ( @(x) islogical(x) || iscell(x) ) );
%
parse(p,varargin{:});

%
if isempty(p.Results.maxn)
    maxNCol = length(obj.trialNames);
else
    maxNCol = p.Results.maxn;
end

if iscell(p.Results.scores)
    scoreInd = ismember(scores, p.Results.scores);
else
    scoreInd = p.Results.scores;
end

%%
copyObj = obj.deepCopy;

%% prefiltering
ind = false(size(copyObj.cell_ID,1),1);
if ~isempty(p.Results.minspikes)
    tmp                                       = cellfun(@(x) cellfun(@(x) length(x),x,'uni',0),copyObj.spikeData.spk_Times,'uni',0);
    nSpikes                                   = cell2mat(horzcat(tmp{:}));
    ind(all(nSpikes < p.Results.minspikes,2)) = true;
end

% apply filter
copyObj.deleteData('cells',ind);

%% Set up results table %
scoreDum  = nan(1,maxNCol);
varList   =   {
                'cellID',      NaN; ...
                'rat',         NaN; ...
                'age',         scoreDum; ...
                'dataset',     cell(1,1); ...
                %'trialInd',    scoreDum; ...
                'envType',     cell(size(scoreDum)); ...
                
                'posData',     cell(size(scoreDum)); ...

                'nSpks',       scoreDum; ...
                'meanRate',    scoreDum; ...
                %'peakRate',    scoreDum; ...

                'rateMap',     cell(size(scoreDum)); ...
                'posMap',      cell(size(scoreDum)); ...
                'dirMap',      cell(size(scoreDum)); ...

                'spkTimes',    cell(size(scoreDum)); ...

                'SI',          scoreDum; ...
                'RV',          scoreDum; ...
                'gridness',    scoreDum; ...
                'gridness_ell',scoreDum; ...
                'intraStab',   scoreDum; ...
                % 'interStab',   nan; ...        % 
                'borderScore', scoreDum; ...

      };
varList = varList';
%
if strcmp(p.Results.rowformat,'cell')
    nCells                    = size(copyObj.cell_ID,1);
else
    nCells = 1;
end
ResT                          = repmat(cell2table(varList(2,:)),nCells,1); % repmat seems to be the only way to pre-allocate table with n rows and variable column format 
ResT.Properties.VariableNames = varList(1,:);

%% assigning index
if ~isempty(p.Results.tindex)
    [dataInd, tabInd] = getIndex({copyObj.trialMetaData.trialType},p.Results.tindex);
else
    [dataInd, tabInd] = deal(1:length(copyObj.trialNames));
end

prmsRate = copyObj.mapParams.rate;

%%
% populate fields

if isfield(copyObj.trialMetaData,'animal')
    ResT.rat               = copyObj.trialMetaData(1).animal .* ones(nCells,1);
elseif isfield(copyObj.trialMetaData,'anNum')
    ResT.rat              = copyObj.trialMetaData(1).anNum .* ones(nCells,1);
end
if isfield(copyObj.trialMetaData,'age')
    ResT.age(:,tabInd)    = [copyObj.trialMetaData(dataInd).age] .* ones(nCells,1);
else
    ResT                  = removevars(ResT,'age');
end
ResT.dataset              = repmat({copyObj.dataSetName},nCells,1);
if ~all(cellfun('isempty',{copyObj.trialMetaData(dataInd).trialType}))
    ResT.envType(:,tabInd) = repmat({copyObj.trialMetaData(dataInd).trialType},nCells,1);
else
    ResT                   = removevars(ResT,'envType');
end

%
switch p.Results.rowformat
    case 'cell'

        ResT = removevars(ResT,'posData');
        %
        if strcmp(copyObj.fileType,'.set')
            ResT.cellID = copyObj.cell_ID(:,1:2);
        else
            ResT.cellID = copyObj.cell_ID(:,1);
        end
        %
        tmp                     = cellfun(@(x) cellfun(@(x) length(x),x,'uni',0),copyObj.spikeData.spk_Times(dataInd),'uni',0);
        nSpikes                 = cell2mat(horzcat(tmp{:}));
        ResT.nSpks(:,tabInd)    = nSpikes;

        % 
        ResT.meanRate(:,tabInd) = bsxfun(@rdivide,nSpikes,[copyObj.trialMetaData(dataInd).duration]);

        %
        ResT.spkTimes(:,tabInd) = horzcat(copyObj.spikeData.spk_Times{dataInd});
        %
        if p.Results.addmaps
            if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                copyObj.addMaps('rate',dataInd);
            end
            ResT.rateMap(:,tabInd) = horzcat(copyObj.maps.rate{dataInd});
            ResT.posMap(:,tabInd)  = horzcat(copyObj.maps.pos{dataInd});

            if any(cellfun('isempty',copyObj.maps.dir(dataInd)))
                copyObj.addMaps('dir',dataInd);
            end
            ResT.dirMap(:,tabInd)  = horzcat(copyObj.maps.dir{dataInd});
        end

        % loop over trials
        if any(scoreInd)
            c = 1;
            for i = dataInd
                if scoreInd(1)
                    if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                        copyObj.addMaps('rate',dataInd);
                    end
                    ResT.SI(:,tabInd(c))           = scanpix.analysis.spatial_info(copyObj.maps.rate{i},copyObj.maps.pos{i});
                end
                if scoreInd(2)
                    if any(cellfun('isempty',copyObj.maps.dir(dataInd)))
                        copyObj.addMaps('dir',dataInd);
                    end
                    ResT.RV(:,tabInd(c))           = cell2mat( cellfun(@(x) scanpix.analysis.rayleighVect(x),copyObj.maps.dir{i},'uni',0) );
                end
                if scoreInd(3)
                    if any(cellfun('isempty',copyObj.maps.sACs(dataInd)))
                        copyObj.addMaps('sac',dataInd);
                    end
                    tmp                            = cell2mat( cellfun(@(x) scanpix.analysis.gridprops(x,'getell',true),copyObj.maps.sACs{i},'uni',0) );
                    ResT.gridness(:,tabInd(c))     = tmp(:,1);
                    ResT.gridness_ell(:,tabInd(c)) = tmp(:,2);
                end
                if scoreInd(4)
                    if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                        copyObj.addMaps('rate',dataInd);
                    end
                    ResT.borderScore(:,tabInd(c))  = cell2mat( cellfun(@(x) scanpix.analysis.getBorderScore(x,copyObj.mapParams.rate.binSizeSpat),copyObj.maps.rate{i},'uni',0) );
                end
                %
                if scoreInd(5)
                    mapSeries                      = scanpix.maps.makeMapTimeSeries(copyObj,[0 copyObj.trialMetaData(i).duration/2],i,'prms',prmsRate);
                    ResT.intraStab(:,tabInd(c))    = scanpix.analysis.spatialCorrelation(mapSeries{1},mapSeries{2});
                end
                c = c + 1;
            end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    case 'dataset'
        if strcmp(copyObj.fileType,'.set')
            ResT.cellID             = {copyObj.cell_ID(:,1:2)};
        else
            ResT.cellID             = {copyObj.cell_ID(:,1)};
        end
        ResT.posData(1,tabInd)      = copyObj.posData.XY(dataInd);
        tmp                         = cellfun(@(x) cellfun(@(x) length(x),x,'uni',0),copyObj.spikeData.spk_Times,'uni',0);
        nSpikes                     = cell2mat(horzcat(tmp{:,dataInd}));
        ResT.nSpks(1,tabInd)        = tmp(dataInd);
        ResT.meanRate(1,tabInd)     = num2cell(bsxfun(@rdivide,nSpikes,[copyObj.trialMetaData(dataInd).duration]),1);

        %
        ResT.spkTimes(1,tabInd)     = copyObj.spikeData.spk_Times(dataInd);
        %
        if p.Results.addmaps
            if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                copyObj.addMaps('rate',dataInd);
            end
            ResT.rateMap(1,tabInd)  = copyObj.maps.rate(dataInd);
            ResT.posMap(1,tabInd)   = copyObj.maps.pos(dataInd);

            if any(cellfun('isempty',copyObj.maps.dir(dataInd)))
                copyObj.addMaps('dir',dataInd);
            end
            ResT.dirMap(1,tabInd)   = copyObj.maps.dir(dataInd);
        else
            ResT = removevars(ResT,{'rateMap','posMap','dirMap'});
        end

        if any(scoreInd)

            %
            if scoreInd(1)
                if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                    copyObj.addMaps('rate',dataInd);
                end
                ResT.SI(1,tabInd)    = cellfun(@(x,y) scanpix.analysis.spatial_info(x,y),copyObj.maps.rate(dataInd),copyObj.maps.pos(dataInd),'uni',0);
            end
            if scoreInd(2)
                if any(cellfun('isempty',copyObj.maps.dir(dataInd)))
                    copyObj.addMaps('dir',dataInd);
                end
                ResT.RV(1,tabInd)    = cellfun(@(y) cellfun(@(x) scanpix.analysis.rayleighVect(x),y),copyObj.maps.dir(dataInd),'uni',0);
            end
            if scoreInd(3)
                if any(cellfun('isempty',copyObj.maps.sACs(dataInd)))
                    copyObj.addMaps('sac',dataInd);
                end
                tmp                         = cellfun(@(y) cell2mat(cellfun(@(x) scanpix.analysis.gridprops(x,'getell',true),y,'uni',0)),copyObj.maps.sACs(dataInd),'uni',0);
                ResT.gridness(1,tabInd)     = cellfun(@(x) x(:,1),tmp,'uni',0);
                ResT.gridness_ell(1,tabInd) = cellfun(@(x) x(:,2),tmp,'uni',0);
            end
            if scoreInd(4)
                if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                    copyObj.addMaps('rate',dataInd);
                end
                ResT.borderScore(1,tabInd)  = cellfun(@(y) cellfun(@(x) scanpix.analysis.getBorderScore(x,copyObj.mapParams.rate.binSizeSpat),y),copyObj.maps.rate(dataInd),'uni',0);
            end
            if scoreInd(5)
                tmpStab                     = cell(1,length(dataInd));
                c = 1;
                for i = dataInd
                    mapSeries              = scanpix.maps.makeMapTimeSeries(copyObj,[0 copyObj.trialMetaData(i).duration/2],i);
                    tmpStab(:,c)           = {scanpix.analysis.spatialCorrelation(mapSeries{1},mapSeries{2})};
                    c = c + 1;
                end
                ResT.intraStab(1,tabInd)   = tmpStab;
            end

        end
end
%
ResT = removevars(ResT,scores(~scoreInd));
if ~ismember('gridness', ResT.Properties.VariableNames)
    ResT = removevars(ResT,'gridness_ell');
end

end

% -------------------------------------------------------------------------------------------------
% --- INLINE FUNCTIONS ----------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
function [dataIndex, tableIndex] = getIndex(trialSeqObj,pattern2extract)  





end