function scalePosition(obj, trialIndex, mode, envSz, minOccForEdge)
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
%       scanpix.maps.scalePosition(obj,trialIndices)
%       scanpix.maps.scalePosition(obj,trialIndices, minOccForEdge)
%       scanpix.maps.scalePosition(obj,trialIndices, minOccForEdge, boxExtent)
%
% Inputs:
%    obj           - dacq or npix class object
%    trialIndices  - numeric index of trials or 'all' 
%    minOccForEdge - min dwell time in s for defining edge. first/last pixels that 
%                    that satisfy this is x and y for min/max extent (default=1)
%    boxExtent     - [n1 n2]; number of pixels in env in x and y after scaling
%
% Outputs:
%
% see also:
%
% LM/TW 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input
if nargin < 2
    uiInput = inputdlg({'trialIndex','dwell thresh','Env. size X','Env. size Y (optional)'},'',1,{'','1','62.5',''});
    if isempty(uiInput)
        warning('scaNpix::Maps::scalePosition:You don''t seem to like scaled positions. Have to abort here...');
        return;
    else
        trialIndex = str2double(uiInput{1});
        minOccForEdge = str2double(uiInput{2});
        envSz = [obj.trialMetaData(trialIndex).ppm / 100 * str2double(uiInput{3}), obj.trialMetaData(trialIndex).ppm / 100 * str2double(uiInput{4})];
    end
end

% assume defaults is 'minOccForEdge' and/or 'boxExtent' aren't supplied
if nargin < 4
    envSz = [250 250];
    if strcmp(mode,'envsampling')
        minOccForEdge = 1;
    end
end

if length(envSz) == 1 || isnan(envSz(2))
    envSz = [envSz(1) envSz(1)];
end

% scale path
XYScaled = nan((size(obj.posData.XY{trialIndex})));  % don't use raw pos??


switch lower(mode)
    case 'envsampling'
        
        minOccForEdge = minOccForEdge * obj.params('posFs'); % samples
        
        for j = 1:2
            
            tempPos = obj.posData.XY{trialIndex}(:,j);
            
            pathHist = histcounts( tempPos, 0.5:1:max(obj.posData.XY{trialIndex}(:))  );    % is this right? 768 is the maximum extent for the DACQ camera.
            
            lowerEdge = find(pathHist >= minOccForEdge, 1, 'first');
            upperEdge = find(pathHist >= minOccForEdge, 1, 'last');
            
            tempPos( tempPos > upperEdge )  = NaN;
            tempPos( tempPos <= lowerEdge ) = NaN;       % Doing <=lowerEdge, then subtracting lowerEdge (line 23), makes the lower limit zero, and therefore the first pixel 1.something.
            tempPos = tempPos - lowerEdge;
            tempPos = tempPos .* (  envSz(j) / (upperEdge-lowerEdge) );
            
            tempPos(tempPos > envSz(j)) = envSz(j);    % Have checked that when this happens, it is a many decimal places rounding error.
            
            XYScaled(:,j) = tempPos;
        end
        
    case 'envdimensions'
        
        tempPos = obj.posData.XY{trialIndex};
        tempPos = tempPos - obj.trialMetaData(trialIndex).envBorderCoords(:,1)';
        tempPos( tempPos <= 0 ) = 1;
        for j = 1:2
            upperEdge = diff(obj.trialMetaData(trialIndex).envBorderCoords(j,:));
            tempPos( tempPos(:,j) > upperEdge, j ) = upperEdge;
            
            tempPos(:,j) = tempPos(:,j) .* (  envSz(j) / (upperEdge-1) );
            
            XYScaled(:,j) = tempPos(:,j);
        end    
end
% For consistency, make sure that all x=nan and all y= nan match up %
XYScaled(any(isnan(XYScaled),2),:) = nan;
% output
obj.posData.XY{trialIndex} = XYScaled;


end


