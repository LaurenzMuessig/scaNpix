function ResT = addShufData2Table(ResT, ResShuf, scores, options)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%% 
arguments
    ResT {mustBeA(ResT,'table')}
    ResShuf {mustBeA(ResShuf,'table')}
    scores (1,:) {mustBeA(scores,'cell')}   = {'gridness','gridness_ell','SI','RV','borderScore','speedScore'};
    options.minSpikes (1,:) {mustBeNumeric} = 100;
    options.ageBins (:,2) {mustBeNumeric}   = [21 100];
    options.pctls (1,:) {mustBeNumeric}     = [90, 95, 97.5, 99];
    options.useMEConly (1,1) {mustBeNumericOrLogical} = true;
end

%%
nTrials             = size(ResT.age,2);
% for pop thresh we want to exclude low rate cells
nSpks               = ResT.nSpks;
nSpks(isnan(nSpks)) = 0;
minSpikesInd        = nSpks < options.minSpikes;
%
% for gridness we will also generate a threshold based on both regular and elliptic gridness
if any(strcmp(scores,'gridness')) && any(strcmp(scores,'gridness_ell'))  
    combinedGridnessThresh = true;
else
    combinedGridnessThresh = false;
end

%% Add thresholds by cell
for i = 1:length(scores)

    if ~any(strcmp(scores{i},ResT.Properties.VariableNames)) && ~any(strcmp(scores{i},ResShuf.Properties.VariableNames))
        warning('scaNpix::shuffle::addShufData2Table: Can''t find ''%s'' as a variable name in either the Res table or the ResShuf table. You need to work on your coding skills!',scores{i});
        continue
    end
    %
    tmp                             = scanpix.shuffle.getShufThresh(ResShuf.(scores{i}),'type','cell','pctls',options.pctls);
    ResT.([scores{i} '_cellPctls']) = tmp(:,1:nTrials);
    % 
    if combinedGridnessThresh 
        allGridness                 = cellfun(@(x,y) [x y],ResShuf.gridness,ResShuf.gridness_ell,'UniformOutput',false);
        tmp                         = scanpix.shuffle.getShufThresh(allGridness,'type','cell','pctls',options.pctls);
        ResT.gridnessAll_cellPctls  = tmp(:,1:nTrials);
    end
end

%%
if combinedGridnessThresh
    allGridnessPop               = allGridness(:,1:nTrials);
    allGridnessPop(minSpikesInd) = {nan(size(allGridnessPop{1}))};
    if options.useMEConly; allGridnessPop(~ResT.inMEC,:) = {nan(size(allGridnessPop{1}))}; end
end

%% add pop thresholds
for i = 1:length(scores)

    if ~any(strcmp(scores{i},ResT.Properties.VariableNames)) && ~any(strcmp(scores{i},ResShuf.Properties.VariableNames))
        warning('scaNpix::shuffle::addShufData2Table: Can''t find ''%s'' as a variable name in either the Res table or the ResShuf table. You need to work on your coding skills!',scores{i});
        continue
    end
    %
    tmpVals               = ResShuf.(scores{i})(:,1:nTrials);
    tmpVals(minSpikesInd) = {nan(size(tmpVals{1}))};
    if options.useMEConly; tmpVals(~ResT.inMEC,:) = {nan(size(tmpVals{1}))}; end
    %
    pctlVals              = nan(height(ResT),length(options.pctls));
    if combinedGridnessThresh; pctlValsAllGridness = nan(height(ResT),length(options.pctls)); end
    %
    for k = 1:size(options.ageBins,1)
        indAge                                   = ResT.age(:,1) >= options.ageBins(k,1) & ResT.age(:,1) <= options.ageBins(k,2);
        pctlVals(indAge,1:length(options.pctls)) = ones(sum(indAge),1) .* scanpix.shuffle.getShufThresh(tmpVals(indAge,:),'type','pop','pctls',options.pctls);
        %
        if combinedGridnessThresh
            pctlValsAllGridness(indAge,1:length(options.pctls)) = ones(sum(indAge),1) .* scanpix.shuffle.getShufThresh(allGridnessPop(indAge,:),'type','pop','pctls',options.pctls);
        end
    end
    ResT.([scores{i} '_popPctls']) = pctlVals;
    %
    if combinedGridnessThresh
        ResT.gridnessAll_popPctls  = pctlValsAllGridness; 
        combinedGridnessThresh     = false;
    end
end
%
ResT.Properties.UserData.ageBinsShufflePop = options.ageBins;
ResT.Properties.UserData.pctlsShuffle      = options.pctls;
end
