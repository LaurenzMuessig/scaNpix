function saveFigAsPDF( varargin )
% Save a figure as a somewhat formatted PDF with some options
% package: scanpix.helpers
%
% Usage:    scanpix.helpers.saveFigAsPDF
%           scanpix.helpers.saveFigAsPDF(fileName)
%           scanpix.helpers.saveFigAsPDF(fileName, directory);
%           scanpix.helpers.saveFigAsPDF(figure_handle,_);
%
% Inputs:   fileName         - optional filename for pdf (no extension); default: 'myFig.pdf'  
%           directory        - optional path to directory; default: cd\  
%           figure_handle    - optional figure handle 
%
%
% LM 2019/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PARSE INPUTS
defaultFName  = 'myFig';
defaultDir    = [cd filesep];

if isempty(varargin) || ~any(isgraphics(varargin{1}))
    hFig = gcf;
else
    hFig = varargin{1};
    varargin = varargin(2:end);
end
% parse
p = inputParser;
addOptional(p,'fileName',defaultFName,@(s)ischar(s) || @(s)isstring(s));
addOptional(p,'dirOut',defaultDir,@(s)ischar(s) || @(s)isstring(s));
parse(p,varargin{:});

% output file
fNameOut = fullfile(p.Results.dirOut,[p.Results.fileName '.pdf']);
% check if file is open and append filename in case - maybe over the top,
% but getting an error in case you have figure open in Adobe is annoying
addN = 1;
if fopen(fNameOut,'a') == -1
    while fopen(fNameOut,'a') == -1
        fNameOut = fullfile(p.Results.dirOut,[p.Results.fileName '_' num2str(addN) '.pdf']);
        addN     = addN + 1;
    end
end
  
%% 
hSlider = findall(hFig,'style','slider'); % check if plot is scrollable
if ~isempty(hSlider)
    % a bit more involved for scrollable plot
    hCanvas     = findall(hFig,'type','uipanel');
    canvasUnits = get(hCanvas,'units'); % keep record (although safe to assume it's pixels)
    set(hCanvas,'position',[0 0 1 1].*get(hCanvas,'position'));
    set(hCanvas,'units','centimeters');
    figSz       = get(hCanvas,'position');
    set(hCanvas,'units',canvasUnits); % reset units
    % temporarily hide sliders
    hSlider     = findall(hFig,'style','slider');
    set(hSlider,'visible','off');
    % set the paper properties
    set(hFig,'PaperPositionMode','auto','PaperUnits','centimeters','PaperType','<custom>','PaperSize',[ figSz(3:4) ]); 
    set(hFig,'PaperPosition',[0 0 1 1].*get(hFig,'PaperPosition')); % need to make sure we start paper at [0,0]
else
    figureUnits = get(hFig,'units'); % keep record
    set(hFig,'units','centimeters');
    figSz = get(hFig,'position');
    set(hFig,'units',figureUnits); % reset figure units to previous type
    set(hFig,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[ figSz(3:4) ]); 
end

%% PRINT TO FILE
print(hFig,'-painters','-loose','-dpdf',fNameOut);

% restore slider
set(hSlider,'visible','on');

% annoyingly need to close file - could do fclose('all'), but better to
% be safe than sorry..
fIDs      = fopen('all');
fNameOpen = arrayfun(@fopen, fIDs, 'UniformOutput', 0); % all open file IDs
idx       = strcmp(fNameOut,fNameOpen);
arrayfun(@fclose, fIDs(idx));

end

