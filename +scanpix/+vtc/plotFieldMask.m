function plotFieldMask(BFMask,PoPrFMask,ax,options)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

arguments
    BFMask {mustBeNumericOrLogical}
    PoPrFMask {mustBeNumericOrLogical}
    ax {ishghandle(ax, 'axes')} = axes;
    options.maskCM (:,3) {mustBeNumeric} = [0.847058824 0.847058824 0.847058824; 0 0.439215686 0.749019608; 1 0.937254902 0; 0.568627451 0.819607843 0.309803922;1 1 1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots (3): Barrier field map.
% Generate the mask image that shows barrier field, the post-barrier field and the overlap.
plotMask              = ones(size(BFMask));  % 1= background visited
plotMask( BFMask==3 ) = 5;  % 5=Unvisited.
% Fill in unvisited bins in env centre before proceeding.
P = regionprops( plotMask==5, 'Area', 'PixelIdxList' );
for itPx = 1:length(P)
    if P(itPx).Area<=6;  plotMask( P(itPx).PixelIdxList ) = 1;  end
end
plotMask( BFMask==1 )    = 2;  % 2=Barrier field
plotMask( PoPrFMask==1 ) = 3;  % 3=Post-probe field
plotMask( BFMask==1 & PoPrFMask==1 ) = 4; % 4=Overlap between barrier and post-probe fields.

% Crop plotmask in the same way as rate map.
% plotMask = mapTrimSF( plotMask, mapSdLen, trmCrds);
% Plot.
image( ax, plotMask );
colormap( ax, options.maskCM );
axis(ax, 'equal', 'tight', 'off');

end