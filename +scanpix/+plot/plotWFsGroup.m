function plotWFsGroup(waveForms,cell_IDs)
% plot all waveforms (mean+/-std) based on a cell array  with waveform 
% samples of different cells. Can help to check if clusters might have to
% be merged. We'll make a figure that is scrollable.
% Note: Can be slow if you try to plot too many cells. 
%
% Usage:    
%           scanpix.plot.plotWFsGroup(waveForms)
%           scanpix.plot.plotWFsGroup(waveForms,cell_IDs)
%
% Inputs:   
%           waveForms - cell array of waveform samples of different single units
%           cell_IDs  - cell array of cell ID strings (optional)      
%
% Outputs:  
%
% LM 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TO DO:
% add plot as individual waveforms + mean?
% dynamically change position of cell ID label - for npix better in bottom
% of plot
% for npix could use actual channel ID in title rather than just numeric
% count

plotSize       = [100 90]; % Pixel
plotSep        = [20 30]; % Pixel
offsetBase     = [40 40];
figWdth        = 550;
baseHghtFigure = 200;


%% params
if nargin < 2
    cell_IDs = strcat({'cell_'},num2str((1:length(waveForms))'));
end

%% plot
screenSz      = get(0,'screensize');
% in case v small screen
if figWdth > 0.7*screenSz(3)
    figWdth = 0.7*screenSz(3);
end
hghtScaleFact = (length(waveForms)-1) * (baseHghtFigure/1.5); % add this to fig width for any additional cell
figHght       = baseHghtFigure + hghtScaleFact;
% height
if figHght > 0.7*screenSz(4)
    figHght = 0.7*screenSz(4);
end
figSize       = [0.1*screenSz(3) 0.1*screenSz(4) figWdth figHght]; % normalised units
hScroll       = scanpix.plot.createScrollPlot(figSize); % open scrollable figure
offsets       = offsetBase;
for i = 1:length(waveForms)
    % mean +/- STD
    meanWF    = squeeze(nanmean(waveForms{i},1));
    stdWFs    = squeeze(nanstd(waveForms{i},[],1));   
    maxVal    = ceil((1.05 * ceil( max(meanWF(:)) + max(stdWFs(:)) ))/10)*10; % for plot axis limits
    minVal    = floor((1.05 * floor( min(meanWF(:)) - max(stdWFs(:)) ))/10)*10; % for plot axis limits
    for j = 1:size(meanWF,2)
        % add axis
        hAx = scanpix.plot.addAxisScrollPlot( hScroll,[offsets plotSize], plotSep );
        axes(hAx);
        % plot waveforms
        scanpix.fxchange.shadedErrorBar([],meanWF(:,j),stdWFs(:,j),{'r-','linewidth',2,});
        if i == 1 && length(waveForms) ~= 1
            set(gca,'xlim',[0 size(meanWF,1)],'xtick',[],'ylim',[minVal maxVal]);
        elseif i == length(waveForms)  
            set(gca,'xlim',[0 size(meanWF,1)],'xtick',[],'xticklabel',{''},'ylim',[minVal maxVal]);
            text(0.1*max(xlim), 1.3*max(ylim),['Channel_' num2str(j)],'FontSize',12,'Interpreter','none'); 
        else
            set(gca,'xlim',[0 size(meanWF,1)],'xtick',[],'xticklabel',{''},'ylim',[minVal maxVal]);
        end
        text(0.45*max(xlim),0.8*max(ylim),cell_IDs{i},'FontSize',12,'Interpreter','none');

        if j == 1
            set(gca,'ytick',[floor(min(meanWF(:))) 0 ceil(max(meanWF(:)))],'yticklabel',[floor(min(meanWF(:))),0,ceil(max(meanWF(:)))]);
        else
            set(gca,'ytick',[floor(min(meanWF(:))) 0 ceil(max(meanWF(:)))],'yticklabel',{''});
        end

        % update offsets    
        offsets(1) = offsets(1) + plotSize(1) + plotSep(1);
    end
    % update offsets  
    offsets(2)     = offsets(2) + plotSize(2) + plotSep(2);
    offsets(1)     = offsetBase(1);
end

end

