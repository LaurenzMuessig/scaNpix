function ResShuf = generateShuffData(objData,mode,options)
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

%%
arguments
    objData {mustBeA(objData,'cell')}
    mode (1,:) {mustBeMember(mode,{'cell','pop'})} = 'cell';
    options.nshuf (1,1) {mustBeNumeric} = 10000;  % this many shuffles - only relevant for population shuffle
    options.minshift (1,1) {mustBeNumeric} = 10;     % in sec
    options.mindur (1,1) {mustBeNumeric} = 600;    %
    options.minspikes {mustBeScalarOrEmpty} = []; 
    options.age (:,2) {mustBeNumeric} = [20 25; 26 40];
    options.scores {mustBeA(options.scores,'cell')} = {'SI','RV','gridness','borderScore'};
    options.trialSelectMode (1,:) {mustBeMember(options.trialSelectMode ,{'all','type','rand','pattern'})}  = 'all';
    options.trialtype (1,:) {mustBeText} = 'fam';
end
%%

warning('off');
matchflag = false;

%%
% figure out the max trials in the overall data that will be shuffled so
% data all has the same format
if strcmp(options.trialSelectMode,'all')
    maxNtrials  = max(cellfun(@(x) length(x.trialNames), objData));
elseif strcmp(options.trialSelectMode,'rand')
    maxNtrials = 1;
elseif strcmp(options.triatrialSelectModelsel,'type')
    allTrials  = cellfun(@(x) {x.trialMetaData.trialType}, objData, 'UniformOutput',0);
    maxNtrials = max(cellfun(@(x) sum(x), cellfun(@(x) ismember(x,options.trialType), allTrials,'UniformOutput',0)));
elseif strcmp(options.trialSelectMode,'pattern')
    maxNtrials = length(options.trialtype);
end

%%
ResShuf = makeTable(0,maxNtrials); % preallocate

age    = cell2mat(cellfun(@(x) x.trialMetaData(1).age, objData, 'UniformOutput',0))';
nCells = cell2mat(cellfun(@(x) length(x.spikeData.spk_Times{1}),objData,'uni',0))';
nStep  = options.mindur - 2*(options.minshift) + 1;
% for population shuffle figure out how many shifts/cell we need to reach n
% of population
if strcmp(mode,'pop')
    nShiftCellsByAge = nan(size(age,1),1);
    for i = 1:size(options.age,1)
        nShiftCellsByAge(age >= options.age(i,1) & age <= options.age(i,2)) = ceil( options.nshuf / sum(nCells(age >= options.age(i,1) & age <= options.age(i,2))));
    end
elseif strcmp(mode,'cell')
    nShiftCellsByAge = NaN; % needed as dummy for parfor execution
end
% start parallel pool
checkParPool = gcp('nocreate');
if isempty(checkParPool)
    numCores = feature('numcores');
    parpool(numCores);
end

% set up vars for par toolbox
trialSelectMode = options.trialSelectMode;
trialtype       = options.trialtype;
minspikes       = options.minspikes;
minshift        = options.minshift;
scores          = options.scores;

% loop over data
parfor j = 1:length(objData)

    copyObj = objData{j}.deepCopy; % make a copy of the current object

    if copyObj.mapParams.rate.showWaitBar
        [copyObj.mapParams.rate.showWaitBar,copyObj.mapParams.dir.showWaitBar] = deal(false);
    end
    % work out n shifts per cell
    if strcmp(mode,'pop')
        nShiftCells = nShiftCellsByAge(j);
    elseif strcmp(mode,'cell')
        nShiftCells = nStep;
    end
    % select trials for the shuffling
    if strcmp(trialSelectMode,'rand')
%         ind = ismember({copyObj.trialMetaData.trialType},options.trialtype);
%         numInd = find(ind);
%         selTrial = numInd(randperm(length(numInd),1));
        selTrial         = randi(length(copyObj.trialNames),1);
        delInd           = true(size(ind));
        delInd(selTrial) = false;
        if ~all(delInd)
            copyObj.deleteData('trials',copyObj.trialNames(delInd));
        end
    elseif strcmp(trialSelectMode,'type')
        ind         = find(ismember({copyObj.trialMetaData.trialType},trialtype));
%         selTrial = ind(min(min(options.trialtype),length(ind)):min(max(options.trialtype),length(ind)));
        delInd      = true(size(ind));
        delInd(ind) = false;
        if ~all(delInd)
            copyObj.deleteData('trials',copyObj.trialNames(delInd));
        end
    elseif strcmp(trialSelectMode,'pattern')
        [dataInd, missTrials] = scanpix.helpers.matchTrialSeq2Pattern({copyObj.trialMetaData.trialType},trialtype,'exact',matchflag);
        missTrials(missTrials > length(copyObj.trialNames)) = [];
        if ~isempty(missTrials)
            copyObj.deleteData('trials',copyObj.trialNames(missTrials));
        end
        copyObj.reorderData(dataInd); 
    end
    
    % remove low spike n cells if desired
    if ~isempty(minspikes)
        nSpikes = cell2mat(cellfun(@(y) cell2mat(cellfun(@(x) length(x),y,'UniformOutput',0)),copyObj.spikeData.spk_Times,'UniformOutput',0));
        cellInd = any(nSpikes < minspikes,2);
        if any(cellInd)
            copyObj.deleteData('cells',cellInd);
        end
    end
    %
    tmpT = makeTable(size(copyObj.cell_ID,1),maxNtrials);
    % loop over trials
    for k = 1:length(copyObj.trialNames)

        shiftsAll = linspace( minshift,copyObj.trialMetaData(k).duration - minshift, nStep)'; % all possible shifts

        % shift 'nShiftCells' times
        for l = 1:nShiftCells
            if strcmp(mode,'pop')
                % pick a random shift for each cell from the dataset in case of population shuffle
                shiftInd  = randi(length(shiftsAll),size(copyObj.cell_ID,1),1);
                shifts    = shiftsAll(shiftInd); % diffent shift for each cell
            elseif strcmp(mode,'cell')
                shifts    = shiftsAll(l); % for cell shuffle just pick linearly from distribution of all shifts
            end
            % shift spikes
            copyObj.spikeData.spk_Times{k} = scanpix.shuffle.randSpikes(copyObj.spikeData.spk_Times{k},copyObj.trialMetaData(k).duration,shifts);
            % make maps
            if strcmp(copyObj.type,'npix') && (~isKey(copyObj.params,'InterpPos2PosFs') || ~copyObj.trialMetaData(k).log.InterpPos2PosFs)
                warning('scaNpix::shuffle::generateShuffData:Shuffling neuropix data without interpolation of positions to a fixed sampling rate is not recommended. You might be waiting a while for your data, mate...');
            end

            if ~all(strcmp(scores,'RV'))
                copyObj.addMaps('rate',k);
            end
            if any(strcmp(scores,'gridness'))
                copyObj.addMaps('sac',k);
            end
            if any(strcmp(scores,'RV'))
                copyObj.addMaps('dir',k);
            end
            
            % compute scores; could probably be a bit more elegant, but it
            % works
            for m = 1:length(scores)
                if strcmp(scores{m},'gridness')
                    tmp = copyObj.getSpatialProps(scores{m}, k);
                    tmpT.(scores{m})(:,k) = num2cell([vertcat(tmpT.(scores{m}){:,k}),max(tmp(:,1),tmp(:,4),'omitnan')],2);
                else
                    tmpT.(scores{m})(:,k) = num2cell([vertcat(tmpT.(scores{m}){:,k}),copyObj.getSpatialProps(scores{m}, k)],2);
                end
            end

            % reset shift
            copyObj.spikeData.spk_Times{k} = scanpix.shuffle.randSpikes(copyObj.spikeData.spk_Times{k},copyObj.trialMetaData(k).duration,shifts,'resetFlag',true);
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