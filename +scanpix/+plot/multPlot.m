function [axArray, hScroll] = multPlot(axArraySz,options)
% Create an arbitray grid of axes - subplot^2... Note that figure will be
% fully scrollable (vertical and horizontal), so you can go to infinity and 
% beyond.  
% package: scanpix.plot
%
%
% Syntax:
%       axArray = scanpix.maps.multPlot(axArraySz)
%       axArray = scanpix.maps.multPlot(axArraySz, Name-Value comma separated list )
%
% Inputs:
%    axArraySz   - [m n] - grid size for plot
%    options     - name-value parameter pairs
%
% Outputs:
%   axArray     - array of generated axes
%
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 

%% parse input
arguments
    axArraySz (1,2) {mustBeNumeric}
    options.plotsize (1,2) {mustBeNumeric} = [75 75];  % pixel;
    options.plotsep (1,2) {mustBeNumeric} = [50 30];  % pixel;
    options.offset (1,2) {mustBeNumeric} = [60 50];  % pixel;
    options.figname (1,:) {mustBeText} = 'multPlot'; 
    options.polflag {mustBeNumericOrLogical} = false(axArraySz);
end

% some sanity checcks should go here

%% plot
% set up figure 
screenSz = get(0,'screensize');

figStart = [0.1*screenSz(3) 0.1*screenSz(4)];

nPlots   = prod(axArraySz);
% set width
figWidth = min([nPlots,axArraySz(2)]) * options.plotsize(1) + min([nPlots,axArraySz(2)]) * options.plotsep(1)+3*options.offset(1);
if figWidth > 0.7*screenSz(3)
    sliderStepX = 1/((figWidth - 0.7*screenSz(3)) / (options.plotsize(1)+options.plotsep(1))); % this should roughly result in a step size of 1 plot
     if sliderStepX > 1
        sliderStepX = 1;
    end
    figWidth = 0.7*screenSz(3);
else
    sliderStepX = 1;
end
% set height
figHeight = ceil(nPlots/axArraySz(2)) * options.plotsize(2) + (ceil(nPlots/axArraySz(2))-1)*options.plotsep(2)+3*options.offset(2);
if figHeight > 0.7*screenSz(4)
    sliderStepY = 1/((figHeight - 0.7*screenSz(4)) / (options.plotsize(2)+options.plotsep(2))); % this should roughly result in a step size of 1 plot
    if sliderStepY > 1
        sliderStepY = 1;
    end
    figHeight = 0.7*screenSz(4); 
else
    sliderStepY = 1;
end

% open plot
offsets              = options.offset;
hScroll              = scanpix.plot.createScrollPlot( [figStart  figWidth figHeight ]  );
hScroll.hFig.Name    = options.figname;
hScroll.hFig.Visible = 'off'; % hiding figure speeds up plotting by a fair amount - not sure necessary here when just plotting empty axes?
hScroll.hFig.UserData.plotArray = axArraySz;
%
plotCount = 1;
nRowPlots = 0;
axArray = cell(axArraySz);
for i = 1:axArraySz(1)
    
    for j = 1:axArraySz(2)
        
        % plot
        axArray{i,j} = scanpix.plot.addAxisScrollPlot( hScroll, [offsets options.plotsize], options.plotsep, options.polflag(i,j) );
        
        % update
        offsets(1) = offsets(1) + options.plotsize(1) + options.plotsep(1);
        nRowPlots  = nRowPlots + 1; %;
        
        if nRowPlots >= axArraySz(2)
            offsets(1) = options.offset(1);
            offsets(2) = offsets(2) + options.plotsize(2) + options.plotsep(2);
            nRowPlots  = 0;
        end
        %
        plotCount = plotCount + 1;
    end
end
% set slider step
set(hScroll.hSldX,'SliderStep',[sliderStepX 1]);
set(hScroll.hSldY,'SliderStep',[sliderStepY 1]);
%
hScroll.hFig.Visible = 'on';

end



