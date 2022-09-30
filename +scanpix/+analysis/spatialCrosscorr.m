function [spatCorr] = spatialCrosscorr(interpBlRm, interpCompRm, varargin)
% 2D cross-correlation of rate maps.
% package: scanpix.analysis
%
%       C = scanpix.analysis.spatialCrosscorr(map1, map2);
%       C = scanpix.analysis.spatialCrosscorr(map1, map2, 'smooth');
%       C = scanpix.analysis.spatialCrosscorr(map1, map2, 'smooth', filtSize, sigma);
%
% For Autocorrelations, just give map1 = map2
%
% Adding 'smooth' argument will smooth the correlogram. Default smooth is
% gaussian window, 5x5 bins, sigma=1.5. You can also specify these.
%
% Written by Caswell Barry
% small edits by LM @ 2022

prms.removeMinOverlap = 1;
prms.smooth = 0;
prms.hSize = 5;
prms.sigma = 1.5;

%% parse input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This is the template code for name-value list OR struct passing of parameters -- %
if ~isempty(varargin)                                                                %
    if ischar(varargin{1})                                                           %
        for ii=1:2:length(varargin);   prms.(varargin{ii}) = varargin{ii+1};   end   %
    elseif isstruct(varargin{1})                                                     %
        s = varargin{1};   f = fieldnames(s);                                        %
        for ii=1:length(f);   prms.(f{ii}) = s.(f{ii});   end                        %
    end                                                                              %
end                                                                                  %
% ---------------------------------------------------------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Make visited bin filters, mark unvis as 0 %
interpBlRm(interpBlRm==-1) = NaN;       % Catch unvis = -1 maps
interpCompRm(interpCompRm==-1) = NaN;   %

onesInterpBlRm=ones(size(interpBlRm));
nanInterpBlRm=isnan(interpBlRm);
onesInterpBlRm(nanInterpBlRm)=0;
interpBlRm(nanInterpBlRm)=0;

onesInterpCompRm=ones(size(interpCompRm));
nanInterpCompRm=isnan(interpCompRm);
onesInterpCompRm(nanInterpCompRm)=0;
interpCompRm(nanInterpCompRm)=0;

warning('off', 'MATLAB:divideByZero');

% 2D Covariance array %
rm1xrm2=filter2(interpBlRm, interpCompRm, 'full'); % Raw convolution <sum(map1.*map2) at each displacement, zero padded>.

nBins=filter2(onesInterpBlRm, onesInterpCompRm, 'full'); % N overlapping bins at each displacement

if prms.removeMinOverlap
    nBins(nBins<20) = 0; % Remove autocorr bins with small overlaps.
end

sumRm1=filter2(interpBlRm, onesInterpCompRm, 'full'); % Sum of rate in overlapping bins, map 1
sumRm2=filter2(onesInterpBlRm, interpCompRm, 'full'); % Sum of rate in overlapping bins, map 2

covar=(rm1xrm2- ((sumRm1.*sumRm2)./nBins) )./nBins; %True covariance

% Convert covariance into correlation (r) %
interpBlRmSq=interpBlRm.^2;
interpCompRmSq=interpCompRm.^2;
nBinsSq=nBins.^2;

sumRm1Sq=filter2(interpBlRmSq, onesInterpCompRm, 'full');
sumRm2Sq=filter2(onesInterpBlRm, interpCompRmSq, 'full');


stdRm1=(sumRm1Sq./nBins - (sumRm1.^2)./nBinsSq).^0.5;
stdRm2=(sumRm2Sq./nBins - (sumRm2.^2)./nBinsSq).^0.5;

spatCorr=covar./(stdRm1.*stdRm2);

% Smoooth AC map %
if prms.smooth
    spatCorr = imfilter(spatCorr, fspecial('gaussian', prms.hSize, prms.sigma));
end

warning('on', 'MATLAB:divideByZero');

