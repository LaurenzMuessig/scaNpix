function randST = randSpikes(spikeTimes,trialDur,shift,options)
% randSpikes - shift spike times by a set amount for some data shuffling
% package: scanpix.shuffle
%
%
% Syntax:
%       randST = scanpix.shuffle.randSpikes(spikeTimes,trialDur,shift)
%       randST = scanpix.shuffle.randSpikes(spikeTimes,trialDur,shift,resetFlag)
%
% Inputs:
%    spikeTimes - cell array of spike times
%    trialDur   - duration of recording trial in sec
%    shift      - shift in sec; either a siblge value or an array with a
%                 shift for each cell
%    resetFlag  - true/false; resetFlag == true will do reverse shift in case you need to undo a shift which is useful as you normally don't want to overwrite the spike times in an ephys object 
%
%
% Outputs: 
%
%    randST     - shifted spike times
%
% LM 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
arguments
    spikeTimes {mustBeA(spikeTimes,'cell')} 
    trialDur (1,1) {mustBeNumeric} 
    shift (1,1) {mustBeNumeric} 
    options.resetFlag (1,1) {mustBeNumericOrLogical} = false;
    options.circShift (1,1) {mustBeNumericOrLogical} = true;
end


%%
if length(shift) ~= 1 && length(shift) ~= length(spikeTimes)
    error('Supply either a fixed shift for all units or one shift value per unit. Can''t work otherwise, can it?')
end

if length(shift) ~= length(spikeTimes)
    shift = repmat(shift(1),length(spikeTimes),1);
end
shift = num2cell(shift,2);

%%
if ~options.resetFlag
    % shift spiketimes
    randST = cellfun(@(x,y) x + y, spikeTimes, shift,'uni',0);
    if options.circShift
    % wrap around spike times
        randST = cellfun(@(x) x - (x > trialDur) .* trialDur, randST,'uni',0);
        % randST = cellfun(@(x) x + (x < 0) .* trialDur, randST,'uni',0);
    else
        randST = cellfun(@(x) x(x > 0 & x < trialDur), randST,'uni',0);
    end
else
    % undo a shift
    randST = cellfun(@(x,y) x - y, spikeTimes, shift,'uni',0);
    randST = cellfun(@(x) x + (x < 0) .* trialDur, randST,'uni',0);
end

end