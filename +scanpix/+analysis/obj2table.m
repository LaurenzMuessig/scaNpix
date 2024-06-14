function ResT = obj2table(obj,varargin)
% obj2table - convert a scanpix.ephys object to a table 
% package: scanpix.analysis
%
% Syntax:  
%    scanpix.analysis.obj2table(obj)
%    scanpix.analysis.obj2table(obj,key-value pairs)
%
% Inputs:
%    obj       - scanpix.ephys object
%    varargin  -  
%
%
% Outputs:
%    ResT      - table output
%
% LM 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
scores = {'SI','RV','gridness','borderScore','intraStab'};

%% Params
rowFormat    = 'cell';
maxNTrial    = [];
minNSpikes   = [];
trialIndex   = [];
addMaps      = true;
addScores    = true(1,length(scores));
addWFprops   = true;

%
p = inputParser;
addOptional(p, 'rowformat', rowFormat,    ( @(x) mustBeMember(x,{'cell','dataset'}) ));
addParameter(p,'maxn',      maxNTrial,    ( @(x) isscalar(x) || isempty(x) ) );
addParameter(p,'minspikes', minNSpikes,   ( @(x) isscalar(x) || isempty(x) ) );
addParameter(p,'tindex',    trialIndex,   ( @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x)) );
addParameter(p,'addmaps',   addMaps,      @islogical );
addParameter(p,'scores',    addScores,    ( @(x) islogical(x) || iscell(x) ) );
addParameter(p,'addWFprops',addWFprops,    @islogical );
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
if strcmp(p.Results.rowformat,'dataset')
    scoreDum  = cell(1,maxNCol);
else
    scoreDum  = nan(1,maxNCol);
end
varList   =   {
                'rat',         NaN; ...
                'age',         scoreDum; ...
                'dataset',     cell(1,1); ...
                'cellID',      NaN; ...
                %'trialInd',    scoreDum; ...
                'envType',     cell(size(scoreDum)); ...
                'nExp',        scoreDum; ...
                
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
                %
                'waveForms',   cell(size(scoreDum)); ...
                'spikeWidth',  scoreDum; ...
                'meanAC',      scoreDum; ...

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
    [dataInd, missTrials] = scanpix.helpers.matchTrialSeq2Pattern({copyObj.trialMetaData.trialType},p.Results.tindex);
    tabInd = 1:length(p.Results.tindex);
    tabInd(missTrials) = [];
    % ind = 1:length(copyObj.trialNames);
    % ind(dataInd) = [];
    % if ~isempty(dataInd)
    %     copyObj.deleteData('trials',copyObj.trialNames(ind));
    % end
else
    [dataInd, tabInd] = deal(1:length(copyObj.trialNames));
end

prmsRate = copyObj.mapParams.rate;

%%
% populate fields
if isfield(copyObj.trialMetaData,'animal')
    ResT.rat               = copyObj.trialMetaData(1).animal .* ones(nCells,1);
elseif isfield(copyObj.trialMetaData,'anNum')
    ResT.rat               = copyObj.trialMetaData(1).anNum .* ones(nCells,1);
end
if isfield(copyObj.trialMetaData,'age')
    ResT.age(:,tabInd)     = [copyObj.trialMetaData(dataInd).age] .* ones(nCells,1);
else
    ResT                   = removevars(ResT,'age');
end
ResT.dataset               = repmat({copyObj.dataSetName},nCells,1);
if ~all(cellfun('isempty',{copyObj.trialMetaData(dataInd).trialType}))
    ResT.envType(:,tabInd) = repmat({copyObj.trialMetaData(dataInd).trialType},nCells,1);
else
    ResT                   = removevars(ResT,'envType');
end
if isfield(copyObj.trialMetaData,'nExp')
    ResT.nExp(:,tabInd)    = [copyObj.trialMetaData(dataInd).nExp] .* ones(nCells,1);
else
    ResT                   = removevars(ResT,'nExp');
end
%
switch p.Results.rowformat
    case 'cell'

        ResT = removevars(ResT,'posData');
        %
        ResT.cellID = copyObj.cell_ID(:,1:2);
        %
        tmp                     = cellfun(@(x) cellfun(@(x) length(x),x,'uni',0),copyObj.spikeData.spk_Times(dataInd),'uni',0);
        ResT.nSpks(:,tabInd)    = cell2mat(horzcat(tmp{:}));
        %
        ResT.spkTimes(:,tabInd) = horzcat(copyObj.spikeData.spk_Times{dataInd});
        %
        if p.Results.addmaps
            if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                copyObj.addMaps('rate',dataInd);
            end
            ResT.rateMap(:,tabInd) = horzcat(copyObj.maps.rate{dataInd});
            if any(cellfun(@(x) length(x) > 1,copyObj.maps.pos(dataInd)))
                copyObj.addMaps('pos',dataInd);
            end
            ResT.posMap(:,tabInd)  = repmat(horzcat(copyObj.maps.pos{dataInd}),nCells,1); % a bit redundant to have same pos map for every cell but better than to store the adaptively smoothed pos maps

            if any(cellfun('isempty',copyObj.maps.dir(dataInd)))
                copyObj.addMaps('dir',dataInd);
            end
            ResT.dirMap(:,tabInd)  = horzcat(copyObj.maps.dir{dataInd});
        end
        % we need this for getting the true mean rate of each cell
        if copyObj.mapParams.rate.speedFilterFlagRMaps
            speedLims = [copyObj.mapParams.rate.speedFilterLimitLow copyObj.mapParams.rate.speedFilterLimitHigh];
        else
            speedLims = 'none';
        end

        % loop over trials
        c = 1;
        for i = dataInd
            % get true mean rate in case of speed filtering
            if isKey(copyObj.params,'InterpPos2PosFs') && copyObj.params('InterpPos2PosFs')
                posFs  = copyObj.trialMetaData(i).log.InterpPosFs;
            else
                posFs  = copyObj.params('posFs');
            end
            ResT.meanRate(:,tabInd(c))  = scanpix.analysis.getMeanRate(copyObj.spikeData.spk_Times{i},posFs,copyObj.posData.speed{i},speedLims); % NEED TO THINK ABOUT WEHETHER THIS IS A GOOD IDEA? 
            %
            if scoreInd(1)
                if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                    copyObj.addMaps('rate',dataInd);
                end
                ResT.SI(:,tabInd(c))           = copyObj.getSpatialProps('SI', i);
            end
            if scoreInd(2)
                if any(cellfun('isempty',copyObj.maps.dir(dataInd)))
                    copyObj.addMaps('dir',dataInd);
                end
                ResT.RV(:,tabInd(c))           = copyObj.getSpatialProps('RV', i);
            end
            if scoreInd(3)
                if any(cellfun('isempty',copyObj.maps.sACs(dataInd)))
                    copyObj.addMaps('sac',dataInd);
                end
                tmpGridProps                   = copyObj.getSpatialProps('gridprops', i);
                ResT.gridness(:,tabInd(c))     = tmpGridProps(:,1);
                ResT.gridness_ell(:,tabInd(c)) = tmpGridProps(:,2);
            end
            if scoreInd(4)
                if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                    copyObj.addMaps('rate',dataInd);
                end
                ResT.borderScore(:,tabInd(c))  = copyObj.getSpatialProps('bs', i);
            end
            %
            if scoreInd(5)
                mapSeries                      = scanpix.maps.makeMapTimeSeries(copyObj,[0 copyObj.trialMetaData(i).duration/2],i,'prms',prmsRate);
                ResT.intraStab(:,tabInd(c))    = scanpix.analysis.spatialCorrelation(mapSeries{1},mapSeries{2});
            end
            c = c + 1;
        end

        if p.Results.addWFprops
            % spikeProps                = scanpix.analysis.getWaveFormProps(copyObj);
            % ResT.waveForms(1,tabInd)  = copyObj.spikeData.spk_waveforms(dataInd);
            % ResT.spikeWidth(1,tabInd) = num2cell(squeeze(spikeProps(:,3,dataInd)),1);
            % ResT.meanAC(1,tabInd)     = num2cell(squeeze(spikeProps(:,4,dataInd)),1);
        else
            ResT                      = removevars(ResT,'waveForms','spikeWidth','meanAC');
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    case 'dataset'
        %
        ResT.cellID             = {copyObj.cell_ID(:,1:2)};
        ResT.posData(1,tabInd)  = copyObj.posData.XY(dataInd);
        tmp                     = cellfun(@(x) cellfun(@(x) length(x),x,'uni',0),copyObj.spikeData.spk_Times,'uni',0);
        nSpikes                 = cell2mat(horzcat(tmp{:,dataInd}));
        ResT.nSpks(1,tabInd)    = tmp(dataInd);
        ResT.meanRate(1,tabInd) = num2cell(bsxfun(@rdivide,nSpikes,[copyObj.trialMetaData(dataInd).duration]),1);

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
                ResT.SI(1,tabInd)    = num2cell(copyObj.getSpatialProps('SI',dataInd),1); 
            end
            if scoreInd(2)
                if any(cellfun('isempty',copyObj.maps.dir(dataInd)))
                    copyObj.addMaps('dir',dataInd);
                end
                ResT.RV(1,tabInd)    = num2cell(copyObj.getSpatialProps('rv',dataInd),1); 
            end
            if scoreInd(3)
                if any(cellfun('isempty',copyObj.maps.sACs(dataInd)))
                    copyObj.addMaps('sac',dataInd);
                end
                % tmp                         = copyObj.getSpatialProps('gridprops',dataInd); 
                % ResT.gridness(1,tabInd)     = cellfun(@(x) x(:,1),tmp,'uni',0);
                % ResT.gridness_ell(1,tabInd) = cellfun(@(x) x(:,2),tmp,'uni',0);
            end
            if scoreInd(4)
                if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                    copyObj.addMaps('rate',dataInd);
                end
                ResT.borderScore(1,tabInd) = num2cell(copyObj.getSpatialProps('bs', dataInd),1);
            end
            if scoreInd(5)
                tmpStab                    = cell(1,length(dataInd));
                c = 1;
                for i = dataInd
                    mapSeries              = scanpix.maps.makeMapTimeSeries(copyObj,[0 copyObj.trialMetaData(i).duration/2],i);
                    tmpStab(:,c)           = {scanpix.analysis.spatialCorrelation(mapSeries{1},mapSeries{2})};
                    c = c + 1;
                end
                ResT.intraStab(1,tabInd)   = tmpStab;
            end

        end
        
        if p.Results.addWFprops
            spikeProps                = scanpix.analysis.getWaveFormProps(copyObj);
            ResT.waveForms(1,tabInd)  = copyObj.spikeData.spk_waveforms(dataInd);
            ResT.spikeWidth(1,tabInd) = num2cell(squeeze(spikeProps(:,3,dataInd)),1);
            ResT.meanAC(1,tabInd)     = num2cell(squeeze(spikeProps(:,4,dataInd)),1);
        else
            ResT                      = removevars(ResT,'waveForms','spikeWidth','meanAC');
        end
end
%
ResT = removevars(ResT,scores(~scoreInd));
if ~ismember('gridness', ResT.Properties.VariableNames)
    ResT = removevars(ResT,'gridness_ell');
end

end
