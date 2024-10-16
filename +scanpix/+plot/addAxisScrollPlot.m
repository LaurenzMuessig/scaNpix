
function hAx = addAxisScrollPlot( hScroll, plotPos, plotSep, polarFlag )
% Add new axes to scrollable plot, making sure to increase size of canvas
% as we go along.
% package: scanpix.plot
%
%  Usage:   hAx = scanpix.plot.addAxisScrollPlot( hScroll, plotPos, plotSep )      
%
%  Inputs:  
%           hScroll   - handle to scroll plot figure 
%           plotPos   - [left bottom width height] position vector for
%                       current axis in pixel;
%           plotSep   - separation of plots in x+y (in pixel)
%           polarFlag - set to true if you want to create polaraxis
%
% see also: 'scanpix.plot.createScrollPlot' - create figure window
%
%  LM 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent firstAxFlag

if nargin < 4
    polarFlag = false;
end
%% parse input 
% separation of individual plots in x+y
if nargin < 3 
    plotSep(1) = 20; % pixel
    plotSep(2) = 20; % pixel
end
%% draw axes
% check if it's first axis - need to keep a memory of this
if isempty(firstAxFlag); firstAxFlag = true; else; firstAxFlag = false; end 
%
if ~polarFlag
    hAx     = axes('Parent',hScroll.hPan,'Units','pixels','Position',[plotPos(1) plotPos(2)  plotPos(3:4)]);
else
    hAx     = polaraxes('Parent',hScroll.hPan,'Units','pixels','Position',[plotPos(1) plotPos(2) plotPos(3:4)]);
end

% disableDefaultInteractivity( hAx ); % global setting for whole figure is more efficient

%increase panel height to fit all figures if necessary
%current position of axis and canvas
posAx       = get(hAx,'Position'); % current axis size/position
p_hPan      = get(hScroll.hPan, 'Position'); % canvas size/position

% horizontal size + slider
if posAx(1) + posAx(3) > 0.95*p_hPan(3) && hScroll.hPan.UserData(2) > (length(hScroll.hPan.Children)+1)/hScroll.hPan.UserData(1)
    if firstAxFlag
        set(hScroll.hPan, 'Position',[p_hPan(1:2) 1.05*+plotPos(3)+plotSep(1) p_hPan(4)]); % increase canvas size
    else
        if hScroll.hPan.UserData(2) > (length(hScroll.hPan.Children)+1)/hScroll.hPan.UserData(1)
            set(hScroll.hPan, 'Position',[p_hPan(1:2) p_hPan(3)+plotPos(3)+plotSep(1) p_hPan(4)]); % increase canvas size
        else
            set(hScroll.hPan, 'Position',[p_hPan(1:2) p_hPan(3)+0.95*plotPos(3) p_hPan(4)]); % increase canvas size
        end
    end
    p_hPan      = get(hScroll.hPan, 'Position'); %
    set(hScroll.hSldX, 'Max',p_hPan(3),'Enable','on'); % increase slider range
end

% vertical size + slider
if posAx(2) + posAx(4) > 0.95*p_hPan(4)
    if firstAxFlag
        set(hScroll.hPan, 'Position',[p_hPan(1:2) p_hPan(3) 1.05*plotPos(4)+plotSep(2)]);
    else
        if hScroll.hPan.UserData(1)*hScroll.hPan.UserData(2)-hScroll.hPan.UserData(1) > length(hScroll.hPan.Children)
            set(hScroll.hPan, 'Position',[p_hPan(1:2) p_hPan(3) p_hPan(4)+plotPos(4)+plotSep(2)]);
        else
            set(hScroll.hPan, 'Position',[p_hPan(1:2) p_hPan(3) p_hPan(4)+0.95*plotPos(4)]);
        end
    end
    p_hPan      = get(hScroll.hPan, 'Position'); %
    set(hScroll.hSldY,'Max',p_hPan(4),'Value',0,'Enable','on'); % increase range
end


% drawnow

end
