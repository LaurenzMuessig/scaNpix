function [gridPropsOut, spatACsOut] = selectBestGridProps(spatACs,varargin)
% selectBestGridProps - select the best grid properties (using gridness) according to different
% selection criteria
% package: scanpix.analysis
%
% Syntax:  
%    scanpix.analysis.selectBestGridProps(spatACs)
%    scanpix.analysis.selectBestGridProps(spatACs,key-value pairs)
%
% Inputs:
%    spatACs      - cell array of spatial autocorrelations of grid cells
%    varargin     -  
%
%
% Outputs:
%    gridPropsOut - grid props for selection
%    spatACsOut   - cell array of spatial autocorrelations for selection
%
% LM 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMS
trialSel   = 'max';
nTrials    = [];
%
p = inputParser;
addParameter(p, 'trialSel', trialSel, @ischar);
addParameter(p, 'nTrials',  nTrials,  (@(x) isempty(x) || isnumeric(x)) );
%
parse(p,varargin{:});
%
if isempty(p.Results.nTrials)
    nTrials = 1:size(spatACs,2);
else
    nTrials = p.Results.nTrials;
end
%
if length(nTrials) > 1 && isempty(strfind(p.Results.trialSel,'max'))
    error('If you supply data for More than 1 trial you need to select one of the ''max'' selection methods');
end

%%
[~,gridPropsStruct] = cellfun(@(x) scanpix.analysis.gridprops(x,'getellgridness',true),spatACs,'uni',0);
%
tmpGridProps = nan(size(spatACs,1),length(nTrials)*7,2);
for i = 1:2
    tmpGridProps(:,:,i)  = cell2mat(cellfun(@(x) [x.gridness(i)' x.wavelength(i)' x.orientation(i) x.offset(i)' x.ellOrient(1,2)' x.ellAbScale(1,3)' x.ellAbScale(1,4)'],gridPropsStruct(:,nTrials),'uni',0));
end
spatACsReg = cellfun(@(x) x.acReg,gridPropsStruct(:,nTrials),'uni',0);

if length(nTrials) == 1

else
end

switch p.Results.trialSel
    case 'maxabs'
        tmpGridProps = {reshape(tmpGridProps,size(tmpGridProps,1),7,[])};
        spatACsOut   = {horzcat(spatACs,spatACsReg)};
    case 'maxall'
        tmpGridProps = {reshape(tmpGridProps(:,:,1),size(tmpGridProps,1),7,size(spatACs,2)) reshape(tmpGridProps(:,:,2),size(tmpGridProps,1),7,size(spatACs,2))};
        spatACsOut   = {spatACs,spatACsReg};
    case 'maxreg'
        tmpGridProps = {reshape(tmpGridProps(:,:,1),size(tmpGridProps,1),7,size(spatACs,2))};
        spatACsOut   = {spatACs};
    case 'maxell'
        tmpGridProps = {reshape(tmpGridProps(:,:,2),size(tmpGridProps,1),7,size(spatACs,2))};
        spatACsOut   = {spatACsReg};
    case 'reg'
        tmpGridProps = {tmpGridProps(:,:,1)};
        spatACsOut   = {spatACs(:,nTrials)};
    case 'ell'
        tmpGridProps = {tmpGridProps(:,:,2)};
        spatACsOut   = {spatACsReg(:,nTrials)};
    case 'besttrial'
        tmpGridProps = {reshape(tmpGridProps(:,:,1),size(tmpGridProps(:,:,1),1),7,[])};
        % spatACsOut   = {spatACs(:,nTrials)};
end
%
[xr,xc]              = ndgrid(1:size(tmpGridProps{1},1),1:size(tmpGridProps{1},2));
selInd               = cell(size(tmpGridProps));
if ~isempty(strfind(p.Results.trialSel,'max'))
    for i = 1:length(tmpGridProps)
        [~, maxInd]      = max(tmpGridProps{i}(:,1,:),[],3,'omitnan');
        selInd{i}        = sub2ind(size(tmpGridProps{i}),xr,xc,repmat(maxInd,1,size(tmpGridProps{i},2)));
        %
        acInd            = sub2ind(size(spatACsOut{i}),1:size(spatACsOut{i},1),maxInd');
        spatACsOut{i}    = spatACsOut{i}(acInd)';

    end
elseif strcmp(p.Results.trialSel,'besttrial')
        [~, maxInd]      = max(squeeze(mean(tmpGridProps{1}(:,1,:),1,'omitnan')),[],'omitnan');
        selInd{1}        = sub2ind(size(tmpGridProps{1}),xr,xc,repmat(maxInd,size(tmpGridProps{1},1),size(tmpGridProps{1},2)));
        spatACsOut{1}    = horzcat(spatACs(:,maxInd),spatACsReg(:,maxInd));
else
    selInd               = sub2ind(size(tmpGridProps),xr,xc);    
end
%
gridPropsOut = cell(size(tmpGridProps));
for i = 1:length(tmpGridProps)
    gridPropsOut{i}         = tmpGridProps{i}(selInd{i});
end
%
if length(gridPropsOut) == 1
    gridPropsOut = gridPropsOut{1};
    spatACsOut   = spatACsOut{1};
end

end