function cornerPoints = findBoxCorners(corner1, L1, corner2, L2)
% findBoxCorners - recover complete set of corner coordinates of a
% rectangular box from 2 input points (we assume these are diagonally oppo-
% sing corners). Note that this works even when box walls and camera window 
% are misaligned.
%
% package: scanpix.helpers
%
%
% Syntax:
%       cornerPoints = scanpix.helpers.findBoxCorners(cornerNW, L1, cornerSE, L2)
%
% Inputs:
%    corner1  - pix coordinates of corner1 (col vect)
%    L1       - side length 1 in pix 
%    corner2  - pix coordinates of corner2 (col vect)
%    L2       - side length 2 in pix 
%
% Outputs:
%    cornerPoints - pix coordinates of all four corners
%
% LM 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Distance between points
distC1C2 = sqrt( (corner2(1)-corner1(1))^2 + (corner2(2)-corner1(2))^2 );

if ~(distC1C2 < L1+L2)
    % Two intersection points don't exist
    error('No such rectangle ever existed - you fucked up!' );
end

%%
% a 
a = (L1^2 - L2^2 + distC1C2^2) / (2*distC1C2);
% b = (L2^2 - L1^2 + distC1C2^2) / (2*distC1C2);

% h
h=sqrt(L1^2 - a^2);

%% centre of box
boxCentre = [ corner1(1) + (a/distC1C2)*(corner2(1)-corner1(1)) , corner1(2) + (a/distC1C2)*(corner2(2)-corner1(2)) ];

%% all corner points
allCorners = [  corner1,...
                [boxCentre(1) - h*(corner2(2)-corner1(2))/distC1C2 ; boxCentre(2) + h*(corner2(1)-corner1(1))/distC1C2],...
                corner2,...
                [boxCentre(1) + h*(corner2(2)-corner1(2))/distC1C2 ; boxCentre(2) - h*(corner2(1)-corner1(1))/distC1C2] ];
% now sort NW, NE, SE, SW
tmpCorners                  = allCorners;
tmpSums                     = sum(allCorners,1);
[~,SEind]                   = max(tmpSums);
[~,NWind]                   = min(tmpSums);
tmpCorners(:,[SEind,NWind]) = 0;
[~,NEind]                   = max(tmpCorners(1,:));
[~,SWind]                   = max(tmpCorners(2,:));
%
cornerPoints                = [allCorners(:,NWind),allCorners(:,NEind),allCorners(:,SEind),allCorners(:,SWind)];
            
end