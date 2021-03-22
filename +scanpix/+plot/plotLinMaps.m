function plotLinMaps( linMap, varargin )
% plotLinMaps: plot a linearised rate map
%
%
% LM 2021

%% parse
defaultMapType = 'rate';
defaultColMap  = 'parula';
defaulthAx     = [];

p = inputParser;
checkAx = @(x) ishghandle(x, 'axes') || @isempty;
addOptional(p,'ax',defaulthAx, checkAx);
addParameter(p,'type',defaultMapType,@ischar);
addParameter(p,'colmap',defaultColMap,@ischar);
parse(p,varargin{:});

if isempty(p.Results.ax)
    hAx = axes;
else
    hAx = p.Results.ax;
end

%% plot
switch lower(p.Results.type)
    case 'rate'
        %% plot
        area(hAx,1:length(linMap),linMap,'edgecolor','k','facecolor','r');
        maxY = max(linMap)*1.1;
        if maxY == 0
            maxY = 1;
        end
        set(hAx,'xTick','','ytick',[0 maxY],'yticklabel',{'0' sprintf('%2.1f',maxY)},'ylim',[0 maxY],'xlim',[1 length(linMap)]);
        %plot arm transitions
        % c = round( [0.25 0.5 0.75].*size(linMaps{1},2) );  % This marks quarters of track (arms)
        % hold on
        % for k=1:length(c)
        %     plot( [1 1].*c(k), get(gca,'ylim') , 'k:' ,'linewidth',2);
        % end
        % hold off
    case 'cellpos'
        imagesc(hAx,linMap);
        eval(['colormap(hAx,' p.Results.colmap ')']);
        set(hAx,'xTick',[0.5 size(linMap,2)+0.5],'xticklabel',[],'ytick',[0.5 size(linMap,1)-0.5],'yticklabel',[0,size(linMap,1)],'ylim',[0.5 size(linMap,1)+0.5],'xlim',[0.5 size(linMap,2)+0.5]);
    otherwise
        error(['scaNpix::plot::plotLinMap:'  p.Results.type ' is not a valid option for plotting linearised rate maps. Try ''rate''  or ''cellpos''.' ]);
end






end



