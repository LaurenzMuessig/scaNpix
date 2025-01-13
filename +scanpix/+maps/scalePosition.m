function scalePosition(obj, trialIndex, varargin)
% scalePosition - Find the edges of vis env, and scale path 
% package: scanpix.maps
% 
% This is done so rate maps have standard sizes and reflect real size of env 
% (e.g. a standard 62.5cm square box at a standard resolution of 400ppm would
% result in exactly 250 pixels across.
%
%
% Syntax:
%       scanpix.maps.scalePosition(obj)
%       scanpix.maps.scalePosition(obj, trialIndex)
%       scanpix.maps.scalePosition(obj, trialIndex, Name-Value comma separated list)
%
% Inputs:
%    obj         - dacq or npix class object
%    trialIndex  - numeric index of trial to be scaled 
%
% Outputs:
%
% see also:
%
% LM/TW 2020
% LM: 2024: added circle scaling from BVC paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


envSzPix             = [250 250];
minOccForEdge        = 50;
circleFlag           = false;
% these are hard coded for now - not sure they need to be dynamically set
minOccForRad         = 0.5;
cenFindIterLim       = 4;
cenFindOffsetThr     = 3;
radiusEstimateMethod = 'largestWellSampledRadius';   %  'fixedPercentileOfDwell';  

p = inputParser;
addParameter(p,'envszpix',envSzPix,@isnumeric);
addParameter(p,'minoccedge',minOccForEdge,@isscalar);
addParameter(p,'circflag',circleFlag,@islogical);
addParameter(p,'minoccrad',minOccForRad,@isscalar);
addParameter(p,'cenlim',cenFindIterLim,@isscalar);
addParameter(p,'cenoffset',cenFindOffsetThr,@isscalar);
addParameter(p,'radest',radiusEstimateMethod,@ischar);
parse(p,varargin{:});


% parse input
if nargin < 2
    uiInput = inputdlg({'trialIndex','dwell thresh','Env. size X','Env. size Y (optional)','circleFlag'},'',1,{'','50','62.5','','1'});
    if isempty(uiInput)
        warning('scaNpix::Maps::scalePosition:You don''t seem to enjoy scaling positions. Have to abort here...');
        return;
    else
        trialIndex    = str2double(uiInput{1});
        minOccForEdge = str2double(uiInput{2});
        envSzPix      = [obj.trialMetaData(trialIndex).ppm / 100 * str2double(uiInput{3}), obj.trialMetaData(trialIndex).ppm / 100 * str2double(uiInput{4})];
        circleFlag    = logical(str2double(uiInput{5}));
    end
else 
    envSzPix      = p.Results.envszpix;
    minOccForEdge = p.Results.minoccedge;
    circleFlag    = p.Results.circflag;
end

if length(envSzPix) == 1 || isnan(envSzPix(2))
    envSzPix = [envSzPix(1) envSzPix(1)];
end

if circleFlag
    [XYScaled, lowerEdge] = scaleCircEnvs(obj,trialIndex,10,envSzPix(1),minOccForRad,cenFindIterLim,cenFindOffsetThr,radiusEstimateMethod);
else
    [XYScaled, lowerEdge] = scaleRectEnvs(obj,trialIndex,minOccForEdge,envSzPix);
end

% For consistency, make sure that all x=nan and all y= nan match up %
XYScaled(any(isnan(XYScaled),2),:) = NaN;
% output
obj.posData.XY{trialIndex} = XYScaled;
%
obj.trialMetaData(trialIndex).PosIsFitToEnv{1,1} = true;
obj.trialMetaData(trialIndex).PosIsFitToEnv{1,2} = lowerEdge;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [XYScaled, lowerEdge] = scaleRectEnvs(obj,trialIndex,minOccForEdge,envSzPix)
% scale path
XYScaled  = nan(size(obj.posData.XY{trialIndex}));  % 
lowerEdge = nan(1,2);

for j = 1:2

    tempPos = obj.posData.XY{trialIndex}(:,j);
    pathHist = histcounts( tempPos, 0.5:1:round(max(obj.posData.XY{trialIndex}(:))*1.1)  );    % using the max(pos) should make it unviversal between npix and dacq.

    if isempty(obj.trialMetaData(trialIndex).envBorderCoords)
        lowerEdge(j) = find(pathHist >= minOccForEdge, 1, 'first');
        upperEdge    = find(pathHist >= minOccForEdge, 1, 'last');
    else
        % in case we have the edges recorded as well we take those into
        % account as well
        lowerEdge(j) = min([obj.trialMetaData(trialIndex).envBorderCoords(j,:),find(pathHist >= minOccForEdge, 1, 'first')]);
        upperEdge    = max([obj.trialMetaData(trialIndex).envBorderCoords(j,:),find(pathHist >= minOccForEdge, 1, 'last')]); %
    end
    %
    tempPos( tempPos > upperEdge )     = NaN;
    tempPos( tempPos <= lowerEdge(j) ) = NaN;       % Doing <=lowerEdge, then subtracting lowerEdge (line 23), makes the lower limit zero, and therefore the first pixel 1.something.
    tempPos                            = tempPos - lowerEdge(j);
    tempPos                            = tempPos .* ( envSzPix(j) / (upperEdge-lowerEdge(j)) );
    %
    tempPos(tempPos > envSzPix(j))     = envSzPix(j);
    %
    XYScaled(:,j)                      = tempPos;

end

end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [XYScaled, lowerEdge] = scaleCircEnvs(obj,trialIndex,minOccForEdge,envSzPix,minOccForRad,cenFindIterLim,cenFindOffsetThr,radiusEstimateMethod)
%

% First, to define the centre, find the edges at the cardinal compass points %
% Iterative procedure: find centre, define that as 0,0, rotate 45degrees, find again and mean that centre with 0,0, rotate 45, and repeat until centre diverges from zero by less that X threshold.
tempPos = obj.posData.XY{trialIndex};
nRot = 0;
tempCen = [100 100];
while nRot<cenFindIterLim && all(abs(tempCen-0)>cenFindOffsetThr)   % The last term is the key one, checking that the temp centre is converging to zero. nRot is a failsafe againt non-convergence.
    % Find the edges in this orientation %
    tempEdges = nan(2,2);
    for i=1:2

        pathBins = ( min(tempPos(:,i),[],1,'omitnan')-0.5) : 1 : (max(tempPos(:,i),[],1,'omitnan')+0.5 );
        pathHist = histcounts( tempPos, pathBins  );

        tempEdges(i,1) = pathBins( find(pathHist >= minOccForEdge, 1, 'first') );  % Here, we just want something like 5 - 10 samples, to remove reflections etc.
        tempEdges(i,2) = pathBins( find(pathHist >= minOccForEdge, 1, 'last')  );
    end
    tempCen = mean(tempEdges,2);
    % Now rotate by 45 degrees %
    xCen   = tempPos(:,1) - tempCen(1);
    yCen   = tempPos(:,2) - tempCen(2);
    [TH,R] = cart2pol(xCen,yCen);
    TH     = TH + (pi/4);
    [tempPos(:,1),tempPos(:,2)]  = pol2cart( TH, R );
    % Bump the iteration counter %
    nRot = nRot + 1;
end
TH = TH - ((pi/4)*nRot);   % From this point on, the function should be working on


%%% Find the extent of the visited path (2 possible methods) %%%
if strcmp(radiusEstimateMethod,'largestWellSampledRadius')
    % (1) Starting with the largest distance, check dwell at all distances until we find a distance is sufficiently well sampled %
    dists = unique( round(R(~isnan(R))) );
    dists = flip(dists);
    n=1;  rLim = [];
    minOccForRadAbs = (minOccForRad/100) * length(R);
    for i=1:length(dists)
        if sum( round(R(~isnan(R))) == dists(n) ) >= minOccForRadAbs   % In this case, you want to try approx minOccForEdge=0.02*trial_dur (assuming out 1 pixel of 100, 0.99^2 = 0.98, therefore remaining area shoudl be approx 0.02)
            rLim = dists(n);
            break
        else
            n = n+1;
        end
    end
    if isempty(rLim); error('No radius exists with occ >= minOccForRad'); end

elseif strcmp(radiusEstimateMethod,'fixedPercentileOfDwell')
    % (2) take a fixed percentile of path data (arranged centre -> edge) as the definition of the edge. %
    % Has the advantage it will always find something approximately OK, but disadvantege that you are guarenteed to throw away X% data every time %
    if minOccForRad<50;   minOccForRad = 100-minOccForRad;   end  % In case the percentile is specified as in the low tail, but we always want the high tail.
    rLim = prctile(R, minOccForRad);
end

% Scale that the largest well-sampled distance to centre is set to 100 pixels from centre %
RScaled     = R .* ((envSzPix/2)/rLim);  % boxExtent is divided by two because it is specified as a diameter, here we are using it to set radius.
[xCen,yCen] = pol2cart(TH,RScaled);
XYScaled    = [xCen + (envSzPix/2), yCen + (envSzPix/2)];
% Set the pixels outside the radius limit to NaN %
XYScaled(R > rLim,:) = NaN;
lowerEdge            = [min(obj.posData.XY{trialIndex}(R <= rLim,1)) min(obj.posData.XY{trialIndex}(R <= rLim,2))];

% If requested, plot for testing purposes %
% if prms.debugPlot
%     figure;  plot( C.x, C.y, 'k-' );
%     axis equal;
%     hold on;
%     [xRad,yRad] = pol2cart( linspace(0, pi*2, 360), ones(1,360).*rLim );
%     plot( xRad, yRad, 'b-' );
% end


end

