function hScroll = createScrollPlot( figSize )
% Create a scrollable plot. Make figure, add canvas and one horizontal and 
% vertical slider. See subfunctions for callbacks for slider and wheel
% scroll
% Note that CTRL+mousewheel will enforce horizontal scrolling in figure 
% window
%
% package: scanpix.plot
%
%  Usage:   scanpix.plot.createScrollPlot 
%           scanpix.plot.createScrollPlot( figSize )
%
%  Inputs:  
%           figSize - position vector of figure in pixel, 
%                     values for 'Position' property of figure
%                     ( [left bottom width height] ) 
%         
% see also: 'scanpix.plot.addAxisScrollPlot' - add axes to canvas
%
% LM 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sliderSzPix   = 20; % slider width in pixels, (20 seems good size)

%%
if isempty(figSize)
    screenSz  = get(0,'screensize');
    figSize   = [0.1*screenSz(3) 0.1*screenSz(4) 0.5*screenSz(3) 0.5*screenSz(4)]; % default
end

%% create figure, panel, and slider
% main figure
% Note: Matlab's annoying axis toolbars create a big delay after plotting figure - removing toolbar altogether for the whole figure is more efficient than disabling the interactivity of every axis and making its toolbar invisible 
hScroll.hFig  = figure('Menubar','figure', 'Resize','off','Units','pixel', 'Position',figSize,'WindowScrollWheelFcn',@wheelScroll,'KeyReleaseFcn',@keyReleaseFcn,'KeyPressFcn',@keyPressFcn,'ToolBar','none','Renderer','painters');
figSizePix    = get(hScroll.hFig,'Position');

% make canvas
hScroll.hPan  = uipanel('Parent',hScroll.hFig,'Units','pixels','Position',[0 sliderSzPix figSizePix(3)-sliderSzPix figSizePix(4)-sliderSzPix],'tag','canvas','BackgroundColor','w');
set(hScroll.hPan,'UserData',[figSizePix(3) figSizePix(4)] );
% make sliders
hScroll.hSldX = uicontrol('Parent',hScroll.hFig,'Style','slider', 'Enable','off','Units','pixel','Position',[0 0 figSizePix(3)-sliderSzPix sliderSzPix],...
                            'Min',0,'Max',figSizePix(3)-sliderSzPix, 'Value',0,'SliderStep', [.15 1],'String','X','Callback',{@moveSlider, hScroll.hPan},'tag','sliderX');

% vertical                                               
hScroll.hSldY = uicontrol('Parent',hScroll.hFig,'Style','slider','Enable','off','Units','pixel','Position',[figSizePix(3)-sliderSzPix 0 sliderSzPix figSizePix(4)],...
                            'Min',0,'Max',figSizePix(4)-sliderSzPix, 'Value',0,'SliderStep', [.15 1],'String','Y','Callback',{@moveSlider, hScroll.hPan},'tag','sliderY');
% control of scrolling type                            
hScroll.hFig.UserData.keyHeld  = '';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CALLBACKS %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wheelScroll( hFig, evnt )
% Callback for wheel scrolling in scrollable plot.
%
%  Usage:   wheelScroll( handle, event )          
%      
%  LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%


% we want to update the panel and the slider position

% current panel position
hPan         = findall(hFig,'Tag','canvas');
p            = get(hPan, 'Position');
orgPanelSz   = get(hPan, 'UserData');

% if ctrl is pressed scroll horizontally
% ctrlPress    = get(gcbf, 'Currentkey');
if strcmp(hFig.UserData.keyHeld ,'ctrl')
    prctMove = orgPanelSz(1) / p(3) / 4; % scroll by this much % of panel size / event (1/4 of org panel size seems good)
    hSlide   = findall(hFig,'Tag','sliderX');
    type     = 'x';
else
    prctMove = orgPanelSz(2) / p(4) / 4; % scroll by this much % of panel size / event (1/4 of org panel size seems good)
    hSlide   = findall(hFig,'Tag','sliderY');
    type     = 'y';
end

if evnt.VerticalScrollCount < 0 % scroll up/right
    % update panel - make sure we stay within canvas size
    switch type
        case 'x'
            newPos = max([-p(3)+orgPanelSz(1) p(1)-prctMove*p(3)]);
            set(hPan,'Position',[newPos p(2:4)]);  
        case 'y'
            newPos = max([-p(4)+orgPanelSz(2) p(2)-prctMove*p(4)]);
            set(hPan,'Position',[p(1) newPos  p(3:4)]);  
    end
elseif evnt.VerticalScrollCount > 0 % scroll down/left
    % update panel - make sure we stay within canvas size
    switch type
        case 'x'
            newPos = min([0 p(1)+prctMove*p(3)]);
            set(hPan,'Position',[newPos  p(2:4)]); 
        case 'y'
            newPos = min([0 p(2)+prctMove*p(4)]);
            set(hPan,'Position',[p(1) newPos  p(3:4)]); 
    end
elseif evnt.VerticalScrollCount == 0 
    % protect against touch pad scroll issue
    return
end

% update slider
switch type
    case 'x'
        relPosCanvas = abs(newPos)/(p(3)-orgPanelSz(1));
        set(hSlide,'Value',max([0 relPosCanvas*(p(3))-20])); % -20 is slider size in pix
    case 'y'
        relPosCanvas = abs(newPos)/(p(4)-orgPanelSz(2));
        set(hSlide,'Value',max([0 relPosCanvas*(p(4))-20])); % -20 is slider size in pix
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function moveSlider( hSlider, ~, hPan )
% Callback for sliders in scrollable plot.
%
%  Usage:   moveSlider( hSlider,~,hPan )          
%  
%  Inputs:  hSlider   - handle to slider that was moved 
%           hPan      - handle to canvas
%         
%  LM 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% grab a few slider props
% type
sliderType = get(hSlider, 'String');
% value
offset     = get(hSlider, 'Value');

% current panel position
p          = get(hPan, 'Position');
orgPanelSz = get(hPan, 'UserData');

% update panel position
switch sliderType
    case 'X'
        % how much has slider moved 
        relSlidPos = offset/p(3);
        newVal     = -relSlidPos * (p(3)-orgPanelSz(1)); % subtract original canvas size so we don't scroll beyond panel size
        set(hPan,'Position',[newVal p(2:4)]);
    case 'Y'
        % how much has slider moved 
        relSlidPos = offset/p(4);
        newVal     = -relSlidPos * (p(4)-orgPanelSz(2));
        set(hPan,'Position',[p(1) newVal p(3:4)]);    
end  

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keyPressFcn(hFig,evt)
    if strcmp(evt.Key,'control')
        hFig.UserData.keyHeld  = 'ctrl';
    end
end
function keyReleaseFcn(hFig,evt)
    if strcmp(evt.Key,'control')
        hFig.UserData.keyHeld  = '';
    end
end
