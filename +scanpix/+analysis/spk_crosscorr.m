function [xc, lags]=spk_crosscorr(stA,stB,xcBin,xcWin,trDur,varargin)
% Spike train cross- or auto-correlogram.
%
%       [corr timeLags]=spk_crosscorr(spikeTrainA,spikeTrainB,corr_bin,corr_window,trial_length)
%       [corr timeLags]=spk_crosscorr(spikeTrainA,'AC',corr_bin,corr_window,trial_length)
%       [corr timeLags]=spk_crosscorr(spikeTrainA,spikeTrainB,corr_bin,corr_window,trial_length,  .. optionName,  'optionValue', ..)
%
% spikeTrainA and spikeTrainB are *column* vectors of the spike times to correlate, in sec. Please don't supply row vectors, this will break function.
%
% If you supply the string 'AC' in place of spikeTrainB, function will know this is an autocorrelogram, and set the central bin to NaN.
%
% corr_bin is the correlogram bin, in sec.
% corr_window is the correlogram window length (i.e. maximum lag), in sec.
% trial_length is the total trial duration, in sec.
%
% Optional parameters can be supplied as a list of name/value pairs, or as a structure with .field=name, value=parameter:
% 
%   'plot', axisHandle               : Pass an axis handle, will plot into that axis.
%   'norm', normalisationMode        : Matlab normalistion option for the function XCORR, default is 'unbiased'.
%   'spikeProb, logical_flag         : Make the units of the cross-corr correspond to spike probability (i.e. set XCORR normalistion to 'none' and multiply by 100). It is up to the user to set a value of corr_bin that is small enough to make this meaningful.
%   'normaliseToShuf, logical_flag   : Will create a popluation of 'shuffled' cross-corrs (based on shifting and wrapping one of the two spike trains), real cross-corr is then re-expressed in terms of shuffle (see next line for details).
%   'normaliseToShufMode, 'mode'     : 'norm' (divide by mean of shuffle) or 'Z' (subtract mean of shuffle and divide by SD of shuffle). Default = 'norm'.
%   'nShuffles', N                   : How many shuffles, for above normalistion? Default=100.
%   'shuffleMinOffset', t            : Minimum time shift for shuffle, in sec. Default=20.


%%% Parse input %%%
% Function defaults %
prms.plot                = [];
prms.norm                = 'unbiased';
prms.spikeProb           = 0;
prms.nShuffles           = 100;
prms.normaliseToShuf     = 0;
prms.normaliseToShufMode = 'norm';  % 'norm' (divide by mean of shuffle) or 'Z' (subtract mean and divide by SD, of shuffle).
prms.shuffleMinOffset    = 20;
prms.useGPU              = 0;
% Parse optional parameters: these can be a string-value CSL, or a structure.
if nargin > 5 && ischar(varargin{1})
    for ii=1:2:length(varargin)
        prms.(varargin{ii}) = varargin{ii+1};
    end
elseif nargin > 5 && isstruct(varargin{1})
    pIn = varargin{1};  f=fieldnames(pIn);
    for ii=1:length(f)
        prms.(f{ii}) = pIn.(f{ii});
    end
end

% Some parameters might have to be changed depending on the input:
% 1) If we want a spike probability plot, or we are doing a shuffled normalisation, the normalistaion *by N overlapping bins* (i.e. dependent on lag) is set to 'none'. For the latter case, because this is h)ow xcorr2 runs, I would have to write my own normalisation code.
if prms.spikeProb || prms.normaliseToShuf
    prms.norm = 'none';  
end
% 2) Check if autocorr.
isAutoCorr = 0;
if ischar(stB) && strcmp(stB,'AC')
    isAutoCorr = 1;
    stB = stA;
end
% 3) check if bad input: both spike trains are [].
if isempty(stA) && isempty(stB)
    lags = -(xcWin/xcBin) : (xcWin/xcBin);    xc = nan(size(lags'));
    return
end

% Convert spike times to spike-time histograms (do this with one call to ACCUMARRAY, more code but saves time) %
spkTrIndAll = [ ceil( stA/xcBin );  ceil( stB/xcBin )  ];  % Convert spike times to xcorr bins, but both in one column vector.
cellInd     = [ ones(size(stA));  ones(size(stB)).*2   ];  % This index is to keep the spikes from each cell separated.
if prms.useGPU
    spkTrIndAll = gpuArray( spkTrIndAll );
    cellInd     = gpuArray( cellInd );
end 
spkTrHist   = accumarray( [spkTrIndAll, cellInd], 1, [trDur/xcBin, 2] );
spkTrHistA  = spkTrHist(:,1);
spkTrHistB  = spkTrHist(:,2);


% XC %
if ~prms.normaliseToShuf
    [xc, lags]=xcorr( spkTrHistA, spkTrHistB, floor(xcWin/xcBin), prms.norm );
else
    [xc, lags]=xcorr(spkTrHistA,spkTrHistB,prms.norm);   %  When doing shuffles, cannot use the 'maxlags' option, it being incompatible with the use of XCORR2 (which doesn't have such a option, and then the lags would be out of sync).
end
if prms.useGPU
    xc   = gather(xc);
    lags = gather(lags);
end
if prms.spikeProb
    xc = xc./length(stB) .* 100;
end
if isAutoCorr
    xc(lags==0) = NaN;
end


% Select only the relevant short window %
winInd = lags>=-(xcWin/xcBin) & lags<=(xcWin/xcBin);
xc     = xc(winInd);
lags   = lags(winInd);


% Produce a population of spike-time shifted x-corrs, so as to convert the real xcorr to 'chance level' (SDs beyond the mean of the shuffled population) %
if prms.normaliseToShuf && ~isempty(stA) && ~isempty(stB)
    % Create an array of time-shifed spike times (nRow = nSpikes, nCol = nShift).
    spkTrBShuf = repmat( stB, 1, prms.nShuffles );  % Set up the array.
    spkTrBShuf = bsxfun(@plus, spkTrBShuf, linspace( prms.shuffleMinOffset, trDur-prms.shuffleMinOffset, prms.nShuffles ));  % Add a unique shift to each column.
    spkTrBShuf( spkTrBShuf>trDur ) = spkTrBShuf( spkTrBShuf>trDur ) - trDur;  % Wrap around spikes now hanging off end of trial.
    % Convert into a series of bin indices, for feeding into ACCUMARRAY %
    spkTrBShufInd = ceil( spkTrBShuf/xcBin );          % Convert to xcorr bins
    spkTrBShufInd( spkTrBShufInd==0 ) = 1;             % Reset bin = 0 to bin = 1.
    spkTrBShufInd = reshape( spkTrBShufInd, [], 1 );   % Reshape to column vector.
    shufRepInd = reshape( repmat( 1:prms.nShuffles, size(spkTrBShuf,1), 1 ),  [], 1 );   % Create an column vector index (for ACCUMARRAY) that says which shuffle each spike belongs to.
    % Make the 2D (col=shift) spike time histogram %
    spkTrHistBShuf = accumarray( [spkTrBShufInd, shufRepInd], 1, [trDur/xcBin, prms.nShuffles] );
    % Run the xcorr of spike train A versus shuffled spike train B
    xcShuf   = xcorr2(spkTrHistA, spkTrHistBShuf);
    xcShuf = xcShuf( winInd, : );
    % Get the mean and SD, and normalise the real xc against these %
    shufMean = nanmean( xcShuf, 2 );
    if strcmp( prms.normaliseToShufMode, 'norm' )
        xc = xc./shufMean;
    elseif strcmp( prms.normaliseToShufMode, 'Z' )
        shufSD   = nanstd(  xcShuf, 0, 2 );
        xc = (xc-shufMean) ./ shufSD;
    end
end

% If requested, plot %
if ~isempty(prms.plot)
    hAxis=prms.plot;
    xVect = (-xcWin:xcBin:xcWin)-(xcBin/2);
    if prms.normaliseToShuf  && strcmp(prms.normaliseToShufMode,'Z')
        % Z-score normalised plots can go negative, so they are done as line plots %
        plot(hAxis,xVect,xc,'k-');
    else
        % This is the 'classic' style %
        shading(hAxis,'FLAT');
        hBar=bar(hAxis,xVect,xc,'histc');
        set(hBar, 'edgecolor', 'none', 'facecolor', [1 1 1].*0.5);
        yLim=get(hAxis,'ylim');   set(hAxis,'ylim',[0 yLim(2)]);  % Make sure lower y-axis limit os 0.
    end
    if ~prms.normaliseToShuf
        set(hAxis,'ytick',[]);  % For a regular (not normalised) XC, after y-axis lower limit is set to 0, turn axis off. (absloute numbers not important).
    end
    set(hAxis,'xlim',[-xcWin xcWin],'xtick',[-xcWin -xcWin/2 0 xcWin/2 xcWin],'tickdir','out'); 
    if prms.spikeProb ||  prms.normaliseToShuf
        set(hAxis,'ytickmode','auto','ylabelmode','auto')
    end
end




