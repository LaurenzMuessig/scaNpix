function ResT = obj2table(obj,rowFormat,options)
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
arguments
    obj (1,:) {mustBeA(obj,'scanpix.ephys')}
    rowFormat (1,:) {mustBeMember(rowFormat,{'cell','dataset'})} = 'cell';
    options.maxn (1,1) {mustBeNumeric} = length(obj.trialNames); 
    options.minspikes (1,1) {mustBeNumeric};
    options.trialPattern (1,:) {mustBeA(options.trialPattern,'cell')} = {};
    options.scores (1,:) {mustBeA(options.scores,'cell')} = {};
    options.mode (1,:)  {mustBeMember(options.mode,{'pattern','exact'})} = 'exact';
    options.bslKey (1,:) {mustBeA(options.bslKey,'cell')} = {'fam'};
    options.ignKey (1,:) {mustBeA(options.ignKey,'cell')} = {'sleep'};
    options.getFlankBSL (1,1) {mustBeNumericOrLogical} = true;
    options.exactflag (1,1) {mustBeNumericOrLogical} = false;
    options.addmaps (1,1) {mustBeNumericOrLogical} = true;
    options.addgridprops (1,1) {mustBeNumericOrLogical} = true;
    options.addwfprops (1,1) {mustBeNumericOrLogical} = false;
    options.addlfp (1,1) {mustBeNumericOrLogical} = false;
end

%%
scores   = {'SI','RV','gridness','borderScore','intraStab'};
scoreInd = ismember(scores, options.scores);
%  
if options.addlfp && strcmp(rowFormat,'dataset')
    addlfp = true;
else
    addlfp = false;
end
%

%%
copyObj = obj.deepCopy;

%% prefiltering
ind = false(size(copyObj.cell_ID,1),1);
if isfield(options,'minspikes')
    tmp                                       = cellfun(@(x) cellfun(@(x) length(x),x,'uni',0),copyObj.spikeData.spk_Times,'uni',0);
    nSpikes                                   = cell2mat(horzcat(tmp{:}));
    ind(all(nSpikes < options.minspikes,2)) = true;
end

% apply filter
copyObj.deleteData('cells',ind);

%% Set up results table %
if strcmp(rowFormat,'dataset')
    scoreDum  = cell(1,options.maxn);
else
    scoreDum  = nan(1,options.maxn);
end
varList   =   {
                'rat',          NaN; ...
                'age',          nan(1,options.maxn); ...
                'dataset',      cell(1,1); ...
                'cellID',       NaN; ...
                %'trialInd',    scoreDum; ...
                'envType',      cell(size(scoreDum)); ...
                'nExp',         nan(1,options.maxn); ...
                'isPreProbe',   scoreDum; ...
                
                'posData',      cell(size(scoreDum)); ...
                'dirData',      cell(size(scoreDum)); ...
                'speed',        cell(size(scoreDum)); ...

                'nSpks',        scoreDum; ...
                'meanRate',     scoreDum; ...
                %'peakRate',    scoreDum; ...

                'rateMap',      cell(size(scoreDum)); ...
                'posMap',       cell(size(scoreDum)); ...
                'dirMap',       cell(size(scoreDum)); ...

                'linMap',       cell(size(scoreDum)); ...
                'linPos',       cell(size(scoreDum)); ...

                'spkTimes',     cell(size(scoreDum)); ...

                'SI',           scoreDum; ...
                'RV',           scoreDum; ...
                'gridness',     scoreDum; ...
                'gridness_ell', scoreDum; ...
                'gridScale',    scoreDum; ...
                'gridOr',       scoreDum; ...
                'gridScale_ell',scoreDum; ...
                'gridOr_ell',   scoreDum; ...
                'intraStab',    scoreDum; ...
                % 'interStab',   nan; ...        % 
                'borderScore',  scoreDum; ...

                %
                'objectPos',    cell(size(scoreDum)); ...
                %
                'waveForms',    cell(size(scoreDum)); ...
                'spikeWidth',   scoreDum; ...
                'meanAC',       scoreDum; ...
                'lfp',          cell(size(scoreDum)); ...
                'lfpChannel',   scoreDum; ...
                'lfpFilter',    scoreDum; ...

      };
varList = varList';
%
if strcmp(rowFormat,'cell')
    nCells                    = size(copyObj.cell_ID,1);
else
    nCells = 1;
end
ResT                          = repmat(cell2table(varList(2,:)),nCells,1); % repmat seems to be the only way to pre-allocate table with n rows and variable column format 
ResT.Properties.VariableNames = varList(1,:);

%% assigning index
if ~isempty(options.trialPattern)
    [dataInd, missTrials] = scanpix.helpers.matchTrialSeq2Pattern({copyObj.trialMetaData.trialType},options.trialPattern,'bslKey',options.bslKey,'exactflag',options.exactflag,'ignKey',options.ignKey,'getFlankBSL', options.getFlankBSL,'mode',options.mode);
    dataInd               = dataInd(1,:);
    tabInd                = 1:length(dataInd)+length(missTrials);
    tabInd(missTrials)    = [];
else
    [dataInd, tabInd]     = deal(1:length(copyObj.trialNames));
end

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

if isfield(copyObj.trialMetaData,'is_preprobe')
    ResT.isPreProbe(:,tabInd)     = [copyObj.trialMetaData(dataInd).is_preprobe] .* ones(nCells,1);
else
    ResT                   = removevars(ResT,'isPreProbe');
end
if isfield(copyObj.trialMetaData,'nExp')
    ResT.nExp(:,tabInd)    = [copyObj.trialMetaData(dataInd).nExp] .* ones(nCells,1);
else
    ResT                   = removevars(ResT,'nExp');
end

%
% TODO: case when there is no high sample rate lfp for dacq data
if addlfp
    if strcmp(copyObj.fileType,'.set')
        ResT.lfp(:,tabInd)        = copyObj.lfpData.lfpHighSamp(dataInd); 
    
        tmpChanInf = cell(1,length(dataInd));
        for i = 1:length(dataInd)
            tmpChanInf{i} = [copyObj.trialMetaData(dataInd(i)).lfp_channel; copyObj.lfpData.lfpTet{dataInd(i)}];
        end
        ResT.lfpChannel(:,tabInd) = tmpChanInf; 
        ResT.lfpFilter(:,tabInd)  = {copyObj.trialMetaData(dataInd).lfp_filter};
    else
        %%%%%% ?????? %%%%%%
    end
else
    ResT = removevars(ResT,{'lfp','lfpChannel','lfpFilter'});
end

%%
switch rowFormat
    case 'cell'

        ResT                    = removevars(ResT,{'posData','dirData','speed'});
        %
        ResT.cellID             = copyObj.cell_ID(:,1:2);
        %
        tmp                     = cellfun(@(x) cellfun(@(x) length(x),x,'uni',0),copyObj.spikeData.spk_Times(dataInd),'uni',0);
        ResT.nSpks(:,tabInd)    = cell2mat(horzcat(tmp{:}));
        %
        ResT.spkTimes(:,tabInd) = horzcat(copyObj.spikeData.spk_Times{dataInd});
        %
        if options.addmaps
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

             % lin maps are special case - need to be pre-made to avoideovercomplicating things
            if any(strcmpi(ResT.envType(1,tabInd),'sqtrack'))
                % NEEDS WORK
                % idx                = find(strcmpi(ResT.envType(1,tabInd),'sqtrack'));
                % ResT.linMap(:,idx) = horzcat(copyObj.maps.lin(dataInd(idx)));
                % ResT.linPos(:,idx) = horzcat(copyObj.maps.linPos(dataInd(idx)));
            else
                ResT = removevars(ResT,{'linMap','linPos'});
            end
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
                tmpGridProps                        = copyObj.getSpatialProps('gridprops', i);
                ResT.gridness(:,tabInd(c))          = tmpGridProps(:,1);
                ResT.gridness_ell(:,tabInd(c))      = tmpGridProps(:,4);
                if options.addgridprops
                    ResT.gridScale(:,tabInd(c))     = tmpGridProps(:,2);
                    ResT.gridOr(:,tabInd(c))        = tmpGridProps(:,3);
                    ResT.gridScale_ell(:,tabInd(c)) = tmpGridProps(:,5);
                    ResT.gridOr_ell(:,tabInd(c))    = tmpGridProps(:,6);
                end
            end
            if scoreInd(4)
                if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                    copyObj.addMaps('rate',dataInd);
                end
                ResT.borderScore(:,tabInd(c))  = copyObj.getSpatialProps('bs', i);
            end
            %
            if scoreInd(5)
                mapSeries                      = scanpix.maps.makeMapTimeSeries(copyObj,[0 copyObj.trialMetaData(i).duration/2],i);
                ResT.intraStab(:,tabInd(c))    = scanpix.analysis.spatialCorrelation(mapSeries{1},mapSeries{2});
            end

            if isfield(copyObj.trialMetaData,'objectPos') && ~isempty(copyObj.trialMetaData(i).objectPos)
                if copyObj.trialMetaData(i).PosIsScaled
                    copyObj.trialMetaData(i).objectPos = copyObj.trialMetaData(i).objectPos .* (copyObj.trialMetaData(i).ppm/copyObj.trialMetaData(i).ppm_org);
                end
                
                if copyObj.trialMetaData(i).PosIsFitToEnv{1,1}
                    tempObjPos = reshape(copyObj.trialMetaData(i).objectPos,[2 size(copyObj.trialMetaData(i).objectPos,2)/2]);
                    tempObjPos = tempObjPos - [copyObj.trialMetaData(i).PosIsFitToEnv{1,2}(1); copyObj.trialMetaData(i).PosIsFitToEnv{1,2}(2)];
                else
                    tempObjPos = copyObj.trialMetaData(i).objectPos';
                end
                ResT.objectPos(:,tabInd(c)) = {tempObjPos(:)'};
            end
            c = c + 1;
        end
        %
        if all(cellfun('isempty',ResT.objectPos(:)))
            ResT                      = removevars(ResT,{'objectPos'});
        end
        %
        if options.addwfprops
            % spikeProps                = scanpix.analysis.getWaveFormProps(copyObj);
            % ResT.waveForms(1,tabInd)  = copyObj.spikeData.spk_waveforms(dataInd);
            % ResT.spikeWidth(1,tabInd) = num2cell(squeeze(spikeProps(:,3,dataInd)),1);
            % ResT.meanAC(1,tabInd)     = num2cell(squeeze(spikeProps(:,4,dataInd)),1);
        else
            ResT                      = removevars(ResT,{'waveForms','spikeWidth','meanAC'});
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    case 'dataset'
        %
        ResT.cellID             = {copyObj.cell_ID(:,1:2)};
        ResT.posData(1,tabInd)  = copyObj.posData.XY(dataInd);
        ResT.dirData(1,tabInd)  = copyObj.posData.direction(dataInd);
        ResT.speed(1,tabInd)    = copyObj.posData.speed(dataInd);
        tmp                     = cellfun(@(x) cellfun(@(x) length(x),x,'uni',0),copyObj.spikeData.spk_Times,'uni',0);
        nSpikes                 = cell2mat(horzcat(tmp{:,dataInd}));
        ResT.nSpks(1,tabInd)    = tmp(dataInd);
        ResT.meanRate(1,tabInd) = num2cell(bsxfun(@rdivide,nSpikes,[copyObj.trialMetaData(dataInd).duration]),1);

        %
        ResT.spkTimes(1,tabInd)     = copyObj.spikeData.spk_Times(dataInd);
        %
        if options.addmaps
            if any(cellfun('isempty',copyObj.maps.rate(dataInd)))
                copyObj.addMaps('rate',dataInd);
            end
            ResT.rateMap(1,tabInd)  = copyObj.maps.rate(dataInd);
            ResT.posMap(1,tabInd)   = copyObj.maps.pos(dataInd);

            if any(cellfun('isempty',copyObj.maps.dir(dataInd)))
                copyObj.addMaps('dir',dataInd);
            end
            ResT.dirMap(1,tabInd)   = copyObj.maps.dir(dataInd);

            % lin maps are special case - need to be pre-made to avoideovercomplicating things
            if any(strcmpi(ResT.envType(1,tabInd),'sqtrack'))
                idxOut                = find(strcmpi(ResT.envType,'sqtrack'));
                idxIn                 = find(strcmpi({copyObj.trialMetaData.trialType},'sqtrack')); 
                ResT.linMap(1,idxOut) = copyObj.maps.lin(idxIn);
                ResT.linPos(1,idxOut) = copyObj.posData.linXY(idxIn);
            else
                ResT               = removevars(ResT,{'linMap','linPos'});
            end
        else
            ResT                   = removevars(ResT,{'rateMap','posMap','dirMap','linMap','linPos'});
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
                tmp                         = copyObj.getSpatialProps('gridprops',dataInd); 
                % ResT.gridness(1,tabInd)     = cellfun(@(x) x(:,1),tmp,'uni',0);
                % ResT.gridness_ell(1,tabInd) = cellfun(@(x) x(:,2),tmp,'uni',0);
                ResT.gridness(1,tabInd)     = num2cell(squeeze(tmp(:,1,:)),1);
                ResT.gridness_ell(1,tabInd) = num2cell(squeeze(tmp(:,4,:)),1);
                if options.addgridprops
                    ResT.gridScale(1,tabInd)     = num2cell(squeeze(tmp(:,2,:)),1);
                    ResT.gridOr(1,tabInd)        = num2cell(squeeze(tmp(:,3,:)),1);
                    ResT.gridScale_ell(1,tabInd) = num2cell(squeeze(tmp(:,5,:)),1);
                    ResT.gridOr_ell(1,tabInd)    = num2cell(squeeze(tmp(:,6,:)),1);
                    
                end
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
        
        if options.addwfprops
            spikeProps                = scanpix.analysis.getWaveFormProps(copyObj);
            ResT.waveForms(1,tabInd)  = copyObj.spikeData.spk_waveforms(dataInd);
            ResT.spikeWidth(1,tabInd) = num2cell(squeeze(spikeProps(:,3,dataInd)),1);
            ResT.meanAC(1,tabInd)     = num2cell(squeeze(spikeProps(:,4,dataInd)),1);
        else
            ResT                      = removevars(ResT,{'waveForms','spikeWidth','meanAC'});
        end
end

%%
if ~options.addgridprops || ~ismember('gridness', ResT.Properties.VariableNames)
    ResT = removevars(ResT,{'gridScale','gridScale_ell','gridOr','gridOr_ell'});
end

ResT = removevars(ResT,scores(~scoreInd));
if ~ismember('gridness', ResT.Properties.VariableNames)
    ResT = removevars(ResT,'gridness_ell');
end

end
