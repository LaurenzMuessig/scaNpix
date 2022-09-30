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
%       scanpix.maps.scalePosition(obj,trialIndices,mode)
%       scanpix.maps.scalePosition(obj,trialIndices, minOccForEdge)
%
% Inputs:
%    obj           - dacq or npix class object
%    trialIndices  - numeric index of trials or 'all' 
%    mode          - mode for scaling 'envSampling' 'envDimensions' 
%    envSzPix      - [sizePixX sizePixY] - size of environment in pix for x and y 
%    minOccForEdge - min dwell time in s for defining edge. first/last pixels  
%                    that satisfy this is x and y for min/max extent (default=1)
%
% Outputs:
%
% see also:
%
% LM/TW 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


envSzPix       = [250 250];
minOccForEdge  = 1;
circleFlag     = false;

p = inputParser;
addParameter(p,'envszpix',envSzPix,@isnumeric);
addParameter(p,'minoccedge',minOccForEdge,@isscalar);
addParameter(p,'circflag',circleFlag,@islogical);
parse(p,varargin{:});


% parse input
if nargin < 2
    uiInput = inputdlg({'trialIndex','dwell thresh','Env. size X','Env. size Y (optional)','circleFlag'},'',1,{'','1','62.5','','1'});
    if isempty(uiInput)
        warning('scaNpix::Maps::scalePosition:You don''t seem to enjoy scaling positions. Have to abort here...');
        return;
    else
        trialIndex = str2double(uiInput{1});
        minOccForEdge = str2double(uiInput{2});
        envSzPix = [obj.trialMetaData(trialIndex).ppm / 100 * str2double(uiInput{3}), obj.trialMetaData(trialIndex).ppm / 100 * str2double(uiInput{4})];
        circleFlag = logical(str2double(uiInput{5}));
    end
else 
    envSzPix      = p.Results.envszpix;
    minOccForEdge = p.Results.minoccedge;
    circleFlag    = p.Results.circflag;
end

if length(envSzPix) == 1 || isnan(envSzPix(2))
    envSzPix = [envSzPix(1) envSzPix(1)];
end

% scale path
XYScaled = nan((size(obj.posData.XY{trialIndex})));  % don't use raw pos??
for j = 1:2
    
    tempPos = obj.posData.XY{trialIndex}(:,j);
    
    if isempty(obj.trialMetaData(trialIndex).envBorderCoords) || circleFlag
        pathHist = histcounts( tempPos, 0.5:1:round(max(obj.posData.XY{trialIndex}(:))*1.1)  );    % using the max(pos) should make it unviversal between npix and dacq.
        
        lowerEdge = find(pathHist >= minOccForEdge, 1, 'first');
        upperEdge = find(pathHist >= minOccForEdge, 1, 'last');
    else
        lowerEdge = min(obj.trialMetaData(trialIndex).envBorderCoords(j,:));
        upperEdge = max(obj.trialMetaData(trialIndex).envBorderCoords(j,:)); %
    end
    
    tempPos( tempPos > upperEdge )  = NaN;
    tempPos( tempPos <= lowerEdge ) = NaN;       % Doing <=lowerEdge, then subtracting lowerEdge (line 23), makes the lower limit zero, and therefore the first pixel 1.something.
    tempPos = tempPos - lowerEdge;
    tempPos = tempPos .* (  envSzPix(j) / (upperEdge-lowerEdge) );
    
    tempPos(tempPos > envSzPix(j)) = envSzPix(j);
    
    XYScaled(:,j) = tempPos;
end
        
% For consistency, make sure that all x=nan and all y= nan match up %
XYScaled(any(isnan(XYScaled),2),:) = NaN;
% output
obj.posData.XY{trialIndex} = XYScaled;


end


