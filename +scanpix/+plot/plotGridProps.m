function plotGridProps(autoCorr,options)

%%
arguments
    autoCorr {mustBeNumeric}
    options.peakDetect (1,:) {mustBeMember(options.peakDetect,{'watershed','zscore'})} = 'watershed';
    options.zScoreThr (1,1) {mustBeNumeric} = 1;
    options.minPeakSz (1,1) {mustBeNumeric} = 4;
    options.axArray (1,3) {mustBeA(options.axArray,'cell')} = scanpix.plot.multPlot([1 3],'plotsize',[150 150],'plotsep',[75 40]);
end


%%
[~, Props]    = scanpix.analysis.gridprops(autoCorr,'peakDetect',options.peakDetect,'zScoreThr',options.zScoreThr,'minPeakSz',options.minPeakSz);
[~, ellProps] = scanpix.analysis.gridprops(autoCorr,true,'peakDetect',options.peakDetect,'zScoreThr',options.zScoreThr,'minPeakSz',options.minPeakSz,'plotEllipse',true,'ax',options.axArray{2});

%%
for i = 1:2
    %
    if i == 1
        sac      = autoCorr;
        tmpProps = Props;
        ax       = options.axArray{1};
    else
        sac      = ellProps.acReg;
        tmpProps = ellProps;
        ax       = options.axArray{3}; 
    end
    %
    if all(isnan(sac)); continue; end

    %
    sizeAC                             = size(sac);
    % a. plot the sac with gridnessmask in jet and background in grey
    [sacForeGround,sacBackGround]      = deal(sac);
    sacForeGround(~tmpProps.gridMask)  = NaN;
    sacBackGround(tmpProps.gridMask)   = NaN;
    % plot background first (grey)
    imagesc(ax,sacBackGround);
    colormap(ax,gray);
    % set(ax,'ydir','normal');
    % plot foreground (jet) - we need to add new axis to plot (make the background for that transparent)
    ax2 = axes(get(ax,'parent'),'units','pixels','position',get(ax,'position'));
    imagesc(ax2,sacForeGround,'alphadata',~isnan(sacForeGround)); clim(ax2,[min(sac(:),[],'omitnan'),1]);
    set(ax2,'color','none','visible','off');%,'ydir','normal');
    colormap(ax2,jet);
    %
    hold(ax,'on'); hold(ax2,'on');

    % b. plot 3 grid axes
    centralPoint = ceil(sizeAC/2);
    for k=1:3
        plot(ax,[centralPoint(1) tmpProps.closestPeaksCoord(k,1)],[centralPoint(2) tmpProps.closestPeaksCoord(k,2)],'k','Linewidth',2);
        plot(ax2,[centralPoint(1) tmpProps.closestPeaksCoord(k,1)],[centralPoint(2) tmpProps.closestPeaksCoord(k,2)],'k','Linewidth',2);
    end

    % c. plot white horizontal/vertical lines as cartesian reference frame
    plot(ax,[0.5 sizeAC(2)],centralPoint([1 1],[2 2]),'--w','Linewidth',2);
    plot(ax,[centralPoint(1,1) centralPoint(1,1)],[0.5 sizeAC(2)],'--w','Linewidth',2);
    plot(ax2,[0.5 sizeAC(2)],centralPoint([1 1],[2 2]),'--w','Linewidth',2);
    plot(ax2,[centralPoint(1,1) centralPoint(1,1)],[0.5 sizeAC(2)],'--w','Linewidth',2);

    % d. plot a red curve to show the grid offset - depending on which wall grid is anchored to, we need to bake in a shift for y
    [minOffset,ind] = min(tmpProps.offsetFull);
    %
    if ind == 1
        if tmpProps.orientationFull(ind) > 0
            shift = 0;
        else
            shift = 2*pi-minOffset;
        end
    else
        if tmpProps.orientationFull(ind) > 0
            shift = pi/2 - minOffset;
        else
            shift = 3/2*pi;
        end
    end
    %
    mag   = tmpProps.wavelength .* 0.6;
    th    = shift:0.05:minOffset+shift;
    [x,y] = pol2cart(th,mag);
    %
    plot(ax, x + centralPoint(2), centralPoint(1)-y, '-r','Linewidth',3);
    plot(ax2,x + centralPoint(2), centralPoint(1)-y, '-r','Linewidth',3);

    % e. add some annotations
    hold(ax,'off'); hold(ax2,'off');
    axis(ax,'off');
    lastBins = (fliplr(size(sac)));
    text(ax,lastBins(2)+1, 0.5,num2str(round(tmpProps.gridness(1,1),3)));
    text(ax,lastBins(2)+1, lastBins(2)*0.75,sprintf('%s\n%s','scale:',num2str(round(tmpProps.wavelength,2))));
    text(ax,lastBins(2)+1, lastBins(1),sprintf('%s\n%s','orientation:',num2str(round(tmpProps.orientation*180/pi,2))));
    if ind ~= 1
        text(ax,lastBins(2)/2-3.5, lastBins(1)+5,sprintf('%s\n%s','offset:',num2str(round(minOffset*180/pi,2))));
    else
        text(ax,lastBins(2)+1, lastBins(2)/2,sprintf('%s\n%s','offset:',num2str(round(minOffset*180/pi,2))));
    end
end
options.axArray{1}.Title.String = 'spatial autocorr';
options.axArray{3}.Title.String = 'spatial autocorr - ellipse corrected';
%
end