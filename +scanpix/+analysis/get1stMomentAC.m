function moment = get1stMomentAC(st,acWin,acBin,trDur)
% get1stMomentAC - calculates first moment (i.e. mean value) of Auto-correlogram (in ms).
% package: scanpix.analysis
%
%  Usage:  moment = sanpix.analysis.get1stMomentAC(spikeTrain,ACWin,ACBin,trialLength)
%
%  Inputs: st     - cell array with spike times in sec
%          acWin  - window for autocorrelation in ms (max lag)
%          acBin  - bin size for autocorrelation in ms
%          trDur  - trial duration in sec
%
%  Output: moment - mean of autocorrelogram in ms
%
% See Csicsvari J, Hirase H, Czurkó A, Mamiya A, Buzsáki G.1999 for reference values in HPC.
%
% LM/TW 

%%
arguments
    st {mustBeNumeric}
    acWin (1,1) {mustBeNumeric}
    acBin (1,1) {mustBeNumeric}
    trDur (1,1) {mustBeNumeric}
end

%%
if length(st) == 1 || isempty(st)
    moment = NaN;
    return;
end    

% Make spike train vector %
spkTrHist       = histcounts(st, 0 : (acBin/1000) : trDur)';
% AC %
[ac, lags]      = xcorr(spkTrHist, spkTrHist, (acWin/acBin), 'unbiased');
% Select only the relevant short window %
ac( lags == 0 ) = NaN;
winInd          = lags > 0 & lags <= ( acWin/acBin );
ac              = ac(winInd);
lags            = lags(winInd);
%get first moment of ac (i.e. mean in ms) 
%mean = sum(x(i-end)*p(x(i-end))). add binsize/2 to center over bins
moment = ( lags + ( acBin/2 ) ) * (ac / sum(ac));