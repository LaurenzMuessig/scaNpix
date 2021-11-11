function axArray = multPlot(axArraySz,varargin)
% Create an arbitray grid of axes - subplot^2... Note that figure will be
% fully scrollable (vertical and horizontal), so you can go big here.  
% package: scanpix.plot
%
%
% Syntax:
%       axArray = scanpix.maps.multPlot(axArraySz)
%       axArray = scanpix.maps.multPlot(axArraySz, Name-Value comma separated list )
%
% Inputs:
%    axArraySz   - [m n] - grid size for plot
%    varargin    - name-value parameter pairs
%
% Outputs:
%   axArray     - array of generated axes
%
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 

%% parse input
defaultPlotSize       = [75 75];  % pixel
defaultPlotSep        = [40 30];  % pixel
defaultOffsetBase     = [60 50];  % pixel
defaultFigName        = 'multPlot';  
% 
p = inputParser;
addParameter(p,'plotsize',defaultPlotSize);
addParameter(p,'plotsep',defaultPlotSep);
addParameter(p,'offsetbase',defaultOffsetBase);
addParameter(p,'figname',defaultFigName,@ischar);

parse(p,varargin{:});

% some sanity checcks should go here

%% plot
% set up figure 
screenSz       = get(0,'screensize');

figStart = [0.1*screenSz(3) 0.1*screenSz(4)];

nPlots   = prod(axArraySz);
% set width
figWidth = min([nPlots,axArraySz(2)]) * p.Results.plotsize(1) + min([nPlots,axArraySz(2)]) * p.Results.plotsep(1)+3*p.Results.offsetbase(1);
if figWidth > 0.7*screenSz(3)
    figWidth = 0.7*screenSz(3);
end
% set height
figHeight = ceil(nPlots/axArraySz(2)) * p.Results.plotsize(2) + (ceil(nPlots/axArraySz(2))-1)*p.Results.plotsep(2)+3*p.Results.offsetbase(2);
if figHeight > 0.7*screenSz(4)
    figHeight = 0.7*screenSz(4);   
end

% open plot
offsets              = p.Results.offsetbase;
hScroll              = scanpix.plot.createScrollPlot( [figStart  figWidth figHeight ]  );
hScroll.hFig.Name    = p.Results.figname;
hScroll.hFig.Visible = 'off'; % hiding figure speeds up plotting by a fair amount

plotCount = 1;

nRowPlots = 0;
axArray = cell(axArraySz);
for i = 1:axArraySz(1)
    
    for j = 1:axArraySz(2)
        
        % plot
        axArray{i,j} = scanpix.plot.addAxisScrollPlot( hScroll, [offsets p.Results.plotsize], p.Results.plotsep );
        
        % update
        offsets(1) = offsets(1) + p.Results.plotsize(1) + 1.5*p.Results.plotsep(1);
        nRowPlots = nRowPlots + 1; %;
        
        if nRowPlots >= axArraySz(2)
            offsets(1) = p.Results.offsetbase(1);
            offsets(2) = offsets(2) + p.Results.plotsize(2) + p.Results.plotsep(2);
            nRowPlots = 0;
        end
        
        plotCount = plotCount + 1;
    end
end

hScroll.hFig.Visible = 'on';

end



