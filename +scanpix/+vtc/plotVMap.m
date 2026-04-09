function plotVMap(vMap,vMapExtMask,ax,options)
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

%%
arguments
    vMap {mustBeNumeric}
    vMapExtMask {mustBeNumericOrLogical} 
    ax  {ishghandle(ax, 'axes')} = axes;
    options.nSteps (1,:) {mustBeNumeric} = 10;
    options.isCirc {mustBeNumericOrLogical} = true;
    options.vMapsNCrc (1,1) {mustBeNumeric} = 4;
end

%%
if min(vMap(:),[],'omitnan')<0;   vMap = vMap - min(vMap(:),[],'omitnan');   end

%%
% Need to bin by hand, to accomodate extra background colour (=could have been in vMap, but not sampled).
im                               = ceil( vMap./max(vMap(:),[],'omitnan') .* options.nSteps);
im(im==0)                        = 1;
im                               = im + 2;  % Shift the actual rate image up two steps, make room for grey and white backgrounds.
im(isnan(vMap))                  = 1;       % Outside the circle of the vector map is white
im(isnan(vMap) & vMapExtMask==0) = 2;       % Inside the circle but not sampled is dark grey.
cMap                             = [1 1 1; 0.3 0.3 0.3; jet(options.nSteps)];

% Plot
image( ax, im );
ax.Colormap = cMap;


% Plot supplmentary bits on map.
hold(ax,'on');
if options.isCirc
    % Add 'cross-hair' lines.
    imC   = size(im,1)/2;
    hL    = plot(ax, [1 1].*imC, ax.YLim, 'k-');   hL.LineWidth = 0.1;
    hL    = plot(ax, ax.XLim, [1 1].*imC, 'k-');   hL.LineWidth = 0.1;
    diagL = imC * sin(pi/4);
    hL    = plot(ax, ([1 -1].*diagL)+imC, ([1 -1].*diagL)+imC, 'k-');   hL.LineWidth = 0.1;
    hL    = plot(ax, ([1 -1].*diagL)+imC, ([-1 1].*diagL)+imC, 'k-');   hL.LineWidth = 0.1;
    % Add radial distance circles.
    maxR  = length(im)/2;
    rCrcs = linspace(maxR/options.vMapsNCrc, maxR, options.vMapsNCrc);
    th    = linspace(0,2*pi,180);
    for itCr=1:options.vMapsNCrc
        [xC,yC] = pol2cart( th, ones(size(th)).*rCrcs(itCr) );
        hL      = plot( ax, xC+(size(im,2)/2), yC+(size(im,1)/2), 'k-' );   hL.LineWidth = 0.1;
    end
else
    % Quadrant lines.
    for itLn = 1:3
        plot( ax, (size(vMap,2)/4).*itLn.*[1 1], ax.YLim, 'k:' );
    end
    % % Plot props.
    % plot( ax, pkSubs(2), pkSubs(1), 'kx' );
    % plot( ax, VP.ConvexHull(:,1), VP.ConvexHull(:,2), 'k-' );
end

hold(ax,'off');
% Format & add title.
axis(ax,'off','equal');
% title(ax, sprintf('D=%2.1f, A=%d', pkDist, round(pkAng) ));


end