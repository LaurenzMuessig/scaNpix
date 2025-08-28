function saveFigAsPDF( figHandle, options )
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

arguments
    figHandle {ishghandle} = gcf;
    options.filename (1,:) {mustBeText} = 'myFig'; 
    options.dir (1,:) {mustBeFolder} = cd; 
end

% output file
fNameOut = fullfile(options.dir,[options.filename '.pdf']);
% check if file is open and append filename in case - maybe over the top,
% but getting an error in case you have figure open in Adobe is annoying
fNameOut = scanpix.helpers.checkSaveFile(fNameOut);
  
%% 
hSlider = findall(figHandle,'style','slider'); % check if plot is scrollable
if  matlab.graphics.internal.mlprintjob.containsUIElements(figHandle) %~isempty(hSlider)
    % a bit more involved for scrollable plot
    hCanvas     = findall(figHandle,'type','uipanel');
    canvasUnits = get(hCanvas,'units'); % keep record (although safe to assume it's pixels)
    set(hCanvas,'position',[0 0 1 1].*get(hCanvas,'position'));
    set(hCanvas,'units','centimeters');
    figSz       = get(hCanvas,'position');
    set(hCanvas,'units',canvasUnits); % reset units
    % temporarily hide sliders
    set(hSlider,'visible','off');
    % set the paper properties
    set(figHandle,'PaperPositionMode','auto','PaperUnits','centimeters','PaperType','<custom>','PaperSize',[ figSz(3:4) ]); 
    set(figHandle,'PaperPosition',[0 0 1 1].*get(figHandle,'PaperPosition')); % need to make sure we start paper at [0,0]
    % exportapp(figHandle.Parent,fNameOut);
    % exportgraphics(figHandle,fNameOut,'ContentType','vector')
else
    figureUnits = get(figHandle,'units'); % keep record
    set(figHandle,'units','centimeters');
    figSz = get(figHandle,'position');
    set(figHandle,'units',figureUnits); % reset figure units to previous type
    set(figHandle,'PaperPositionMode','auto','PaperUnits','centimeters','PaperSize',[ figSz(3:4) ]); 
end

%% PRINT TO FILE
print(figHandle,'-vector','-loose','-dpdf',fNameOut);

% restore slider
set(hSlider,'visible','on');

% annoyingly need to close file - could do fclose('all'), but better to
% be safe than sorry..
fIDs      = fopen('all');
fNameOpen = arrayfun(@fopen, fIDs, 'UniformOutput', 0); % all open file IDs
idx       = strcmp(fNameOut,fNameOpen);
arrayfun(@fclose, fIDs(idx));

end

