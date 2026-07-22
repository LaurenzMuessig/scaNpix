function plotLinMaps( linMap, options )
% plotLinMaps: plot a linearised rate map
%
%
% LM 2021

%%
arguments
    linMap 
    options.ax {ishghandle(options.ax, 'axes')} = axes;
    options.type (1,:) {mustBeMember(options.type,{'rate','cellpos'})} = 'rate';
    options.colmap (1,:) {mustBeText} = 'parula';
end

%% plot
switch options.type
    case 'rate'
        %% plot
        area(options.ax,1:length(linMap),linMap,'edgecolor','k','facecolor','r');
        maxY = max(linMap)*1.1;
        if maxY == 0
            maxY = 1;
        end
        set(options.ax,'xTick','','ytick',[0 maxY],'yticklabel',{'0' sprintf('%2.1f',maxY)},'ylim',[0 maxY],'xlim',[1 length(linMap)]);
        %plot arm transitions
        % c = round( [0.25 0.5 0.75].*size(linMaps{1},2) );  % This marks quarters of track (arms)
        % hold on
        % for k=1:length(c)
        %     plot( [1 1].*c(k), get(gca,'ylim') , 'k:' ,'linewidth',2);
        % end
        % hold off
    case 'cellpos'
        linMapsNormed = scanpix.maps.normLinMaps(linMap);
        imagesc(options.ax,linMapsNormed);
        eval(['colormap(hAx,' options.colmap ')']);
        set(options.ax,'xTick',[0.5 size(linMapsNormed,2)+0.5],'xticklabel',[],'ytick',[0.5 size(linMapsNormed,1)-0.5],'yticklabel',[0,size(linMapsNormed,1)],'ylim',[0.5 size(linMapsNormed,1)+0.5],'xlim',[0.5 size(linMapsNormed,2)+0.5]);
    otherwise
        error(['scaNpix::plot::plotLinMap:'  options.type ' is not a valid option for plotting linearised rate maps. Try ''rate''  or ''cellpos''.' ]);
end






end



