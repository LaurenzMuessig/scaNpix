function [smoothedRate,smoothedSpk,smoothedPos,medFiltRad] = adaptiveSmooth(pos, spk, alpha, varargin)
% adaptiveSmooth - Adaptive smoothing of rate maps.
% package: scanpix.maps
%
%       [smoothedRate,smoothedSpk,smoothedPos]=scanpix.maps.adaptiveSmooth(posMap,spkMap,alpha)
%
% Each bin is smoothed using a flat, circular kernal. The circle radius 
% is set for each bin, indivdually, such that 
%
%   radius => alpha ./ ( sqrt(nSpike) .* nDwell )
%
% where nSpike and nDwell are the number of spikes, and the amount of dwell time (in s) within the kernel.
%
% smoothedRate, smoothedSpk, smoothedPos are the smoothed maps (spike and pos maps are smoothed 
% with the same kernal group as for the rate map.
%
% TW

% Check for empty spk maps %
if sum(sum(spk))==0
    smoothedPos=pos;    smoothedPos(pos==0)=nan;
    smoothedSpk=spk;    smoothedSpk(pos==0)=nan;
    smoothedRate=spk;   smoothedRate(pos==0)=nan;
    medFiltRad = nan;
    return
end
if ~isempty(varargin) && strcmp(varargin{1},'newKernel');   legacyKernel = 0;   else    legacyKernel = 0;    end
% Pre-assign output %
smoothedPos=zeros(size(pos));
smoothedSpk=zeros(size(pos));
% Visited env template: use this to get numbers of visited bins in filter at edge of environemnt %
vis=zeros(size(pos));
vis(pos>0)=1;
% Pre-assign map which records which bins have passed %
smoothedCheck=false(size(pos));
smoothedCheck(pos==0)=true; % Disregard unvisited - mark as already done.
% Pre-assign list of radii used (this is for reporting purposes, not used for making maps) %
radiiUsedList=nan(1,sum(sum(pos>0)));
radiiUsedCount=1;
% These parameters depend on place or dir mode %
if size(pos,2)>1
    boundary=0;             % IMFILTER boundary condition
    if legacyKernel   
        rBump=0.5;              % Increase radius in 0.5 bin steps.
    else
        rBump=1;
    end
elseif size(pos,2)==1
    boundary='circular';
    rBump=1;                % Increase radius in 1 bin steps.
end


%%% Run increasing radius iterations %%%
r=1; % Circle radius
while any(any(~smoothedCheck))
    % Check radius isn't getting too big (if >map/2, stop running) %
    if r>max(size(pos))/2
        smoothedSpk(~smoothedCheck) = 0;  % Comment LM: This seems the wrong behaviour as NaN is reserved for unvisited positions. Should be 0
        smoothedPos(~smoothedCheck) = eps; % this can't be 0
        break
    end
    % Construct filter kernel ...
    if size(pos,2)>1
        % Place: Flat disk, where r>=distance to bin centre %
        if legacyKernel
            f=fspecial('disk',r); 
            f(f>=(max(max(f))/3))=1;
            f(f~=1)=0;
        else
            [xx,yy] = meshgrid( -r:r );
            [~,rho] = cart2pol( xx, yy );
            f       = double(rho<=r);
        end
    elseif size(pos,2)==1 
        % Direction: boxcar window, r bins from centre symmetrically %
        f=ones(1+(r*2),1);
    end     
    % Filter maps (get N spikes and pos sum within kernel) %
    fSpk=imfilter(spk,f,boundary);
    fPos=imfilter(pos,f,boundary);
    fVis=imfilter(vis,f,boundary);
    % Which bins pass criteria at this radius? %
%     warning('off', 'MATLAB:divideByZero');
    binsPassed=alpha./(sqrt(fSpk).*fPos) <= r;
%     warning('on', 'MATLAB:divideByZero');
    binsPassed=binsPassed & ~smoothedCheck; % Only get the bins that have passed in this iteration.
    % Add these to list of radii used %
    nBins=sum(binsPassed(:));
    radiiUsedList(radiiUsedCount:radiiUsedCount+nBins-1)=r;
    radiiUsedCount=radiiUsedCount+nBins;
    % Assign values to smoothed maps %
    smoothedSpk(binsPassed)=fSpk(binsPassed)./fVis(binsPassed);
    smoothedPos(binsPassed)=fPos(binsPassed)./fVis(binsPassed);
    % Record which bins were smoothed this iteration %
    smoothedCheck(binsPassed)=true;
    % Increase circle radius %
    r=r+rBump;
end

% Assign Output %
% warning('off', 'MATLAB:divideByZero');
smoothedRate=smoothedSpk./smoothedPos;
% warning('on', 'MATLAB:divideByZero');
smoothedRate(pos==0)=nan;
smoothedPos(pos==0) =nan;
smoothedSpk(pos==0) =nan;
medFiltRad = nanmedian(radiiUsedList);

% Report radii sizes %
if 0
    hAllFigs = get(0, 'children');
    hFig = findobj(hAllFigs, 'flat', 'tag', 'adaptiveSmoothPlotWindow');
    if isempty(hFig)
        hFig=figure;
        set(hFig,'tag','adaptiveSmoothPlotWindow');
    else
        figure(hFig);
    end
    hist(radiiUsedList,1:10);
    uiwait(hFig,1.5);
end





