function ResShuf = generateShuffData(objData,varargin)
% generateShuffData - genereate a population of shuffled data score to
% define cell types. Can either do by cell or by age bin (population
% shuffle)
% package: scanpix.shuffle
%
%
% Syntax:
%       ResShuf = scanpix.shuffle.generateShuffData(objData)
%       ResShuf = scanpix.shuffle.generateShuffData(objData,mode)
%       ResShuf = scanpix.shuffle.generateShuffData(__,Name-Value comma separated list)
%
% Inputs:
%    objData     - cell array of scanpix.ephys objects or single object
%    mode        - (optional) - 'cell' or 'pop'
%
%
% Outputs: 
%
%    ResShuf     - Table with shuffled scores
%
% LM 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off');

matchflag = false;

%%
mode            = 'cell';
nShuffles       = 10000;  % this many shuffles - only relevant for population shuffle
minShift        = 10;     % in sec
minTrialDur     = 600;    % 
minNspikes      = [];
ageBins         = [20 25; 26 40];
scores          = {'SI','RV','gridness','borderScore'};
trialSelectMode = 'all';
trialType       = 'fam';

%
p = inputParser;
addOptional(p, 'mode'     ,mode,           ( @(x) mustBeMember(x,{'cell','pop'}) ) );
addParameter(p,'nshuf'    ,nShuffles,      @isscalar);
addParameter(p,'minshift' ,minShift,       @isscalar);
addParameter(p,'mindur'   ,minTrialDur,    @isscalar);
addParameter(p,'minspikes',minNspikes,     @isscalar);
addParameter(p,'age'      ,ageBins                  );
addParameter(p,'scores'   ,scores,         @iscell  );
addParameter(p,'trialsel' ,trialSelectMode,(@(x) mustBeMember(x,{'all','type','rand','pattern'}) ) );
addParameter(p,'trialtype',trialType,      (@(x) ischar(x) || iscell(x) ) );

parse(p,varargin{:});
%
if ~iscell(objData); objData = {objData}; end

% figure out the max trials in the overall data that will be shuffled so
% data all has the same format
if strcmp(p.Results.trialsel,'all')
    maxNtrials  = max(cellfun(@(x) length(x.trialNames), objData));
elseif strcmp(p.Results.trialsel,'rand')
    maxNtrials = 1;
elseif strcmp(p.Results.trialsel,'type')
    allTrials  = cellfun(@(x) {x.trialMetaData.trialType}, objData, 'UniformOutput',0);
    maxNtrials = max(cellfun(@(x) sum(x), cellfun(@(x) ismember(x,p.Results.trialType), allTrials,'UniformOutput',0)));
elseif strcmp(p.Results.trialsel,'pattern')
    maxNtrials = length(p.Results.trialtype);
end

%%
ResShuf = makeTable(0,maxNtrials); % preallocate

age    = cell2mat(cellfun(@(x) x.trialMetaData(1).age, objData, 'UniformOutput',0))';
nCells = cell2mat(cellfun(@(x) length(x.spikeData.spk_Times{1}),objData,'uni',0))';
nStep  = p.Results.mindur - 2*(p.Results.minshift) + 1;
% for population shuffle figure out how many shifts/cell we need to reach n
% of population
if strcmp(p.Results.mode,'pop')
    nShiftCellsByAge = nan(size(p.Results.age,1),1);
    for i = 1:size(p.Results.age,1)
        nShiftCellsByAge(i) = ceil( p.Results.nshuf / sum(nCells(age >= p.Results.age(i,1) & age <= p.Results.age(i,2))));
    end
elseif strcmp(p.Results.mode,'cell')
    nShiftCellsByAge = NaN; % needed as dummy for parfor execution
end
% start parallel pool
checkParPool = gcp('nocreate');
if isempty(checkParPool)
    numCores = feature('numcores');
    parpool(numCores);
end
% loop over data
parfor j = 1:length(objData)

    copyObj = objData{j}.deepCopy; % make a copy of the current object

    if copyObj.mapParams.rate.showWaitBar
        [copyObj.mapParams.rate.showWaitBar,copyObj.mapParams.dir.showWaitBar] = deal(false);
    end
    % work out n shifts per cell
    if strcmp(p.Results.mode,'pop')
        nShiftCells = nShiftCellsByAge(age(j) >= p.Results.age(:,1) & age(j) <= p.Results.age(:,2));
    elseif strcmp(p.Results.mode,'cell')
        nShiftCells = nStep;
    end
    % select trials for the shuffling
    if strcmp(p.Results.trialsel,'rand')
%         ind = ismember({copyObj.trialMetaData.trialType},p.Results.trialtype);
%         numInd = find(ind);
%         selTrial = numInd(randperm(length(numInd),1));
        selTrial         = randi(length(copyObj.trialNames),1);
        delInd           = true(size(ind));
        delInd(selTrial) = false;
        if ~all(delInd)
            copyObj.deleteData('trials',copyObj.trialNames(delInd));
        end
    elseif strcmp(p.Results.trialsel,'type')
        ind         = find(ismember({copyObj.trialMetaData.trialType},p.Results.trialtype));
%         selTrial = ind(min(min(p.Results.trialtype),length(ind)):min(max(p.Results.trialtype),length(ind)));
        delInd      = true(size(ind));
        delInd(ind) = false;
        if ~all(delInd)
            copyObj.deleteData('trials',copyObj.trialNames(delInd));
        end
    elseif strcmp(p.Results.trialsel,'pattern')
        [dataInd, missTrials] = scanpix.helpers.matchTrialSeq2Pattern({copyObj.trialMetaData.trialType},p.Results.trialtype,'exact',matchflag);
        missTrials(missTrials > length(copyObj.trialNames)) = [];
        if ~isempty(missTrials)
            copyObj.deleteData('trials',copyObj.trialNames(missTrials));
        end
        copyObj.reorderData(dataInd); 
    end
    
    % remove low spike n cells if desired
    if ~isempty(p.Results.minspikes)
        nSpikes = cell2mat(cellfun(@(y) cell2mat(cellfun(@(x) length(x),y,'UniformOutput',0)),copyObj.spikeData.spk_Times,'UniformOutput',0));
        cellInd = any(nSpikes < p.Results.minspikes,2);
        if any(cellInd)
            copyObj.deleteData('cells',cellInd);
        end
    end
    %
    tmpT = makeTable(size(copyObj.cell_ID,1),maxNtrials);
    % loop over trials
    for k = 1:length(copyObj.trialNames)

        shiftsAll = linspace( p.Results.minshift,copyObj.trialMetaData(k).duration - p.Results.minshift, nStep)'; % all possible shifts

        % shift 'nShiftCells' times
        for l = 1:nShiftCells
            if strcmp(p.Results.mode,'pop')
                % pick a random shift for each cell from the dataset in case of population shuffle
                shiftInd  = randi(length(shiftsAll),size(copyObj.cell_ID,1),1);
                shifts    = shiftsAll(shiftInd); % diffent shift for each cell
            else
                shifts    = shiftsAll(l); % for cell shuffle just pick linearly from distribution of all shifts
            end
            % shift spikes
            copyObj.spikeData.spk_Times{k} = scanpix.shuffle.randSpikes(copyObj.spikeData.spk_Times{k},copyObj.trialMetaData(k).duration,shifts);
            % make maps
            if strcmp(copyObj.fileType,'.ap.bin') && (~isKey(copyObj.params,'InterpPos2PosFs') || ~copyObj.trialMetaData(k).log.InterpPos2PosFs)
                warning('scaNpix::shuffle::generateShuffData:Shuffling neuropix data without interpolation of positions to a fixed sampling rate is not recommended. You might be waiting a while for your data, mate...');
            end

            if ~all(strcmp(p.Results.scores,'RV'))
                copyObj.addMaps('rate',k);
            end
            if any(strcmp(p.Results.scores,'gridness'))
                copyObj.addMaps('sac',k);
            end
            if any(strcmp(p.Results.scores,'RV'))
                copyObj.addMaps('dir',k);
            end
            
            % compute scores; could probably be a bit more elegant, but it
            % works
            for m = 1:length(p.Results.scores)
                if strcmp(p.Results.scores{m},'gridness')
                    tmp = copyObj.getSpatialProps(p.Results.scores{m}, k);
                    tmpT.(p.Results.scores{m})(:,k) = num2cell([vertcat(tmpT.(p.Results.scores{m}){:,k}),max(tmp(:,1),tmp(:,4),'omitnan')],2);
                else
                    tmpT.(p.Results.scores{m})(:,k) = num2cell([vertcat(tmpT.(p.Results.scores{m}){:,k}),copyObj.getSpatialProps(p.Results.scores{m}, k)],2);
                end
            end

            % reset shift
            copyObj.spikeData.spk_Times{k} = scanpix.shuffle.randSpikes(copyObj.spikeData.spk_Times{k},copyObj.trialMetaData(k).duration,shifts,true);
        end
    end
    % fill in rest of table
    tmpT.rat     = repmat(copyObj.trialMetaData(1).animal,height(tmpT),1);
    tmpT.dataset = repmat({copyObj.dataSetName},height(tmpT),1);
    tmpT.age     = repmat(copyObj.trialMetaData(1).age,height(tmpT),1);
    tmpT.cellID  = copyObj.cell_ID(:,1);
    % final output
    ResShuf      = vertcat(ResShuf,tmpT);
    %
    fprintf('Finished shuffling %s\n',copyObj.dataSetName);
end

% save
save(['ResShuf_' char(datetime('today','format','yyMMdd')) '.mat'],'ResShuf','-v7.3');
%
warning('on');

end

% -------------------------------------------------------------------------------------------------
% --- INLINE FUNCTIONS ----------------------------------------------------------------------------
% -------------------------------------------------------------------------------------------------
function T = makeTable(nCells,maxTrials)

scoreDum  = nan(1,maxTrials);
varList   =   {
    'rat',         NaN; ...
    'dataset',     cell(1,1); ...
    'age',         NaN; ...
    'cellID',      NaN; ...

    'SI',          cell(size(scoreDum)); ...
    'RV',          cell(size(scoreDum)); ...
    'gridness',    cell(size(scoreDum)); ...
    'borderScore', cell(size(scoreDum)); ...

    };
varList = varList';
%
T                          = repmat(cell2table(varList(2,:)),nCells,1); % repmat seems to be the only way to pre-allocate table with n rows and variable column format
T.Properties.VariableNames = varList(1,:);

end