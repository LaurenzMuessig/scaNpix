function cornerPoints = findBoxCorners(cornerNW, L1, cornerSE, L2)


%% Distance between points
distC1C2 = sqrt( (cornerSE(1)-cornerNW(1))^2 + (cornerSE(2)-cornerNW(2))^2 );

if ~(distC1C2 < L1+L2)
    % Two intersection points don't exist
    error('No such rectangle ever existed - you fucked up!' );
end

%% a and b
a = (L1^2 - L2^2 + distC1C2^2) / (2*distC1C2);
% b = (radius2^2 - radius1^2 + distP1P2^2) / (2*distP1P2);

%% h
h=sqrt(L1^2 - a^2);

%% P5
boxCentre = [ cornerNW(1) + (a/distC1C2)*(cornerSE(1)-cornerNW(1)) , cornerNW(2) + (a/distC1C2)*(cornerSE(2)-cornerNW(2)) ];

%% Intersection points
cornerPoints = [ cornerNW,...
                [boxCentre(1) - h*(cornerSE(2)-cornerNW(2))/distC1C2 ; boxCentre(2) + h*(cornerSE(1)-cornerNW(1))/distC1C2],...
                cornerSE,...
                [boxCentre(1) + h*(cornerSE(2)-cornerNW(2))/distC1C2 ; boxCentre(2) - h*(cornerSE(1)-cornerNW(1))/distC1C2] ];
end
