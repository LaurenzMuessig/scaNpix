function cornerPoints = findBoxCorners(cornerNW, L1, cornerSE, L2)
% findBoxCorners - recover complete set of corner coordinates of a
% rectangular box from 2 input points (we assume these are NW and SE
% corners). Note that this works even when box walls and camera window are
% misaligned.
%
% package: scanpix.helpers
%
%
% Syntax:
%       cornerPoints = scanpix.helpers.findBoxCorners(cornerNW, L1, cornerSE, L2)
%
% Inputs:
%    cornerNW - pix coordinates of NorthWest corner (col vect)
%    L1       - side length 1 in pix 
%    cornerSE - pix coordinates of SouthEast corner (col vect)
%    L2       - side length 2 in pix 
%
% Outputs:
%    cornerPoints - pix coordinates of all four corners
%
% LM 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Distance between points
distC1C2 = sqrt( (cornerSE(1)-cornerNW(1))^2 + (cornerSE(2)-cornerNW(2))^2 );

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
boxCentre = [ cornerNW(1) + (a/distC1C2)*(cornerSE(1)-cornerNW(1)) , cornerNW(2) + (a/distC1C2)*(cornerSE(2)-cornerNW(2)) ];

%% all corner points
cornerPoints = [ cornerNW,...
                [boxCentre(1) - h*(cornerSE(2)-cornerNW(2))/distC1C2 ; boxCentre(2) + h*(cornerSE(1)-cornerNW(1))/distC1C2],...
                cornerSE,...
                [boxCentre(1) + h*(cornerSE(2)-cornerNW(2))/distC1C2 ; boxCentre(2) - h*(cornerSE(1)-cornerNW(1))/distC1C2] ];
end