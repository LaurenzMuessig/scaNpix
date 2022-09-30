function xRecFiltInd = findCrossRecCells(obj, trialInd, maps, varargin)
% findCrossRecCells - find duplicate cells that were recorded across tetrodes.
% Method for dacq class objects
% We've been doing this since the Muessig et. al (2019) paper as typically up to 10% of units might be duplicates. 
% We will remove units when pairs have a high spatial correlation and also have a big peak at 0 lag in their Xcorrelogram.
% Note that at the moment you need to explicitly supply maps, we don't check in function whether maps exist in object (this could be changed but
% comes at cost of some complications)
% Also note that by default we will only generate an index instead of removing
% the data from the object
%
% Syntax:  obj.findCrossRecCells(trialInd,maps)
%          obj.findCrossRecCells(trialInd,maps, optionalInputStruct );
%          obj.findCrossRecCells(trialInd,maps, 'inputName', inputVal, .. etc .. );
%
% Inputs:
%    trialInd - index of trial to base check on - this should be a RUN trial and not e.g. a sleep trial as we will use the spatial similarity across cell pairs
%    maps     - optional cell array of rate maps. If not supllied we will make them on the fly. Note that linear rate maps currently need to be supplied as a cell array 
%               for each map (so you need to do 'num2cell(dacqObj.linMaps.linRate{index}, 2)') 
%    varargin - optional prmsStruct or Name-Value pair list to change default values for params 
%
% Outputs:
%
%
% See also:
%
% TW/LM 2019-20


%% PARAMS
prms.xcZeroBinThr        = 3;   % this many more counts in central bin compared to average of rest of cross-correlogram (default=3)
prms.spatialCorrThresh   = 0.8; % min spatial correlation for a cell pair to be considered a duplicate (default=0.8)
prms.ccWin               = 50;  % window size crosscorrelograms in ms (default=50ms)
prms.ccBin               = 1;   % bin size crosscorrelograms in ms (default=1ms)
prms.removeCellsFromData = 0;   % flag if cells should be dropped from data (default=false)

% only relevant for making rate maps on the fly
prms.mapParams           = scanpix.maps.defaultParamsRateMaps; % supply a struct of map params if you want maps to be generated with custom parameters


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for i=1:2:length(varargin);   prms.(varargin{i}) = varargin{i+1};   end      %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for i=1:length(f);   prms.(f{i}) = s.(f{i});   end                           %  
    end                                                                              %
end                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% pre-processing
% make rate maps if not supplied
if nargin < 3 || isempty(maps)
    if isempty(regexp(lower(obj.trialMetaData(trialInd).trialType),'track', 'once') )
        scanpix.maps.addMaps(obj,'rate', trialInd, prms.mapParams.rate);
        maps = obj.maps.rate{trialInd}; % we want each map as a cell array so it's compatible with 2D map format
    else
        scanpix.maps.addMaps(obj,'lin', trialInd, prms.mapParams.lin);
        maps = num2cell(obj.maps.lin{trialInd}{1}, 2); % we want each map as a cell array so it's compatible with 2D map format
    end
end

%% The real thing
STarray = obj.spikeData.spk_Times{trialInd};

% get the relevant cell pair that we want to test, based on their spatial
% similarity
r                = scanpix.analysis.spatialCorrelation(maps);
[rowInd, colInd] = find(triu(r) >= prms.spatialCorrThresh ); % only du upper triangular portion 
diagInd          = rowInd == colInd; % remove diagonals (can't just do r == 1 as there are sometimes many decimal rounding errors)
rowInd           = rowInd(~diagInd);
colInd           = colInd(~diagInd);

% Calculate 0-bin x-correlogram of cells, not on same tetrode, looking for cross-recorded cells %
xcZeroBin = nan( length(rowInd), 1 );
cellPairID = nan( length(rowInd), 2 );
for i = 1:length(rowInd)
    xcTemp                              = scanpix.analysis.spk_crosscorr( STarray{rowInd(i)}, STarray{colInd(i)}, prms.ccBin/1000, prms.ccWin/1000, obj.trialMetaData(trialInd).duration );
    centreInd                           = false(1, length(xcTemp));
    centreInd((length(xcTemp)-1)/2 + 1) = true;
    xcZeroBin(i,1)                      = xcTemp(centreInd) ./ mean( xcTemp( ~centreInd ) );
    cellPairID(i,:)                     = [rowInd(i) colInd(i)];   
end

% Remove cells defined as cross-tetrode recorded - this is slightly involved, as we need to 
% be careful we don't remove more cells than necessary, so go down the list and remove those 
% cells that appear in the most pairs.
xRecPair = cellPairID( xcZeroBin > prms.xcZeroBinThr, : );
cellToRemList = [];
for i = 1:size(xRecPair,1)
   if isnan(xRecPair(i, 1));   continue;   end
   if sum( xRecPair(:) == xRecPair(i,1) )  >=  sum( xRecPair(:) == xRecPair(i, 2) )
       cellToRem = xRecPair(i, 1);
   else
       cellToRem = xRecPair(i, 2); 
   end   
   xRecPair( any(xRecPair == cellToRem,2), : ) = nan;
   cellToRemList = [cellToRemList, cellToRem];
end
% final index
xRecFiltInd                  = true( size( STarray ) );
xRecFiltInd( cellToRemList ) = false;
% obj.uniqueCellIndex          = xRecFiltInd;

% remove cells from object if desired
if prms.removeCellsFromData
    scanpix.helpers.deleteData(obj,'cells',~xRecFiltInd);
end

end

