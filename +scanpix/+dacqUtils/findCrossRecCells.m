function xRecFiltInd = findCrossRecCells(obj, trialInd, options)
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

%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialInd (1,1) {mustBeNumeric}
    options.mapType (1,:) {mustBeMember(options.mapType,{'rate','lin'})} = 'rate';
    options.xcZeroBinThr (1,1) {mustBeNumeric} = 3;
    options.spatialCorrThresh (1,1) {mustBeNumeric} = 0.8; 
    options.ccWin (1,1) {mustBeNumeric} = 50; % in degrees
    options.ccBin (1,1) {mustBeNumeric} = 1; % in degrees
    options.removeCellsFromData (1,1) {mustBeNumericOrLogical} = false; 
end

%% pre-processing
% make rate maps if not supplied
if isempty(obj.maps.(options.mapType){trialInd})
    obj.addMaps(options.mapType, trialInd, obj.mapParams);
    %
    if strcmp(options.mapType,'rate')
        maps = obj.maps.rate{trialInd};
    else
        maps = num2cell(obj.maps.lin{trialInd}{1}, 2);
    end
end

%% The real thing
STarray = obj.spikeData.spk_Times{trialInd};

% get the relevant cell pair that we want to test, based on their spatial
% similarity
r                = scanpix.analysis.spatialCorrelation(maps);
[rowInd, colInd] = find(triu(r) >= options.spatialCorrThresh ); % only du upper triangular portion 
diagInd          = rowInd == colInd; % remove diagonals (can't just do r == 1 as there are sometimes many decimal rounding errors)
rowInd           = rowInd(~diagInd);
colInd           = colInd(~diagInd);

% Calculate 0-bin x-correlogram of cells, not on same tetrode, looking for cross-recorded cells %
xcZeroBin  = nan( length(rowInd), 1 );
cellPairID = nan( length(rowInd), 2 );
for i = 1:length(rowInd)
    xcTemp                              = scanpix.analysis.spk_crosscorr( STarray{rowInd(i)}, STarray{colInd(i)}, options.ccBin/1000, options.ccWin/1000, obj.trialMetaData(trialInd).duration );
    centreInd                           = false(1, length(xcTemp));
    centreInd((length(xcTemp)-1)/2 + 1) = true;
    xcZeroBin(i,1)                      = xcTemp(centreInd) ./ mean( xcTemp( ~centreInd ) );
    cellPairID(i,:)                     = [rowInd(i) colInd(i)];   
end

% Remove cells defined as cross-tetrode recorded - this is slightly involved, as we need to 
% be careful we don't remove more cells than necessary, so go down the list and remove those 
% cells that appear in the most pairs.
xRecPair = cellPairID( xcZeroBin > options.xcZeroBinThr, : );
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
if options.removeCellsFromData
    scanpix.helpers.deleteData(obj,'cells',~xRecFiltInd);
end

end

