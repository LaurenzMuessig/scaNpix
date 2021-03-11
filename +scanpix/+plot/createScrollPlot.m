function hScroll = createScrollPlot( figSize )
% Create a scrollable plot. Make figure, add canvas and one horizontal and 
% vertical slider. See subfunctions for callbacks for slider and wheel
% scroll
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
hScroll.hFig  = figure('Menubar','figure', 'Resize','off','Units','pixel', 'Position',figSize,'WindowScrollWheelFcn',@wheelScroll);
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
                            
% we want record of original canvas height for mouse scrolling limit
set(hScroll.hFig,'UserData',figSizePix(4)); 
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

prctMove = 0.05; % scroll by this much % of panel size / event (5% seems good)

% we want to update the panel and the slider position

% current panel position
hPan        = findall(hFig,'Tag','canvas');
p           = get(hPan, 'Position');
orgPanelSz  = get(hPan, 'UserData');

if evnt.VerticalScrollCount < 0 % scroll up
    % update panel - make sure we stay within canvas size
    newPos  = max([-p(4)+orgPanelSz(2) p(2)-prctMove*p(4)]);
    set(hPan,'Position',[p(1) newPos  p(3:4)]);  
elseif evnt.VerticalScrollCount > 0 % scroll down
    % update panel - make sure we stay within canvas size
    newPos = min([0 p(2)+prctMove*p(4)]);
    set(hPan,'Position',[p(1) newPos  p(3:4)]);  
elseif evnt.VerticalScrollCount == 0 
    % protect against touch pad scroll issue
    return
end

% current slider properties
hSlidY       = findall(hFig,'Tag','sliderY');
% update slider
relPosCanvas = abs(newPos)/(p(4)-orgPanelSz(2));
set(hSlidY,'Value',max([0 relPosCanvas*(p(4))-20])); % -20 is slider size in pix


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

% how much has slider moved 
relSlidPos = offset/p(4);
% update panel position
switch sliderType
    case 'X'
        newVal = -relSlidPos * (p(3)-orgPanelSz(1)); % subtract original canvas size so we don't scroll beyond panel size
        set(hPan,'Position',[newVal p(2:4)]);
    case 'Y'
        newVal = -relSlidPos * (p(4)-orgPanelSz(2));
        set(hPan,'Position',[p(1) newVal p(3:4)]);    
end  

end


