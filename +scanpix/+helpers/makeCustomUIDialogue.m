function output = makeCustomUIDialogue(prompts, defaultVals,fillBoxLength)
% makeCustomUIDialogue - create a custom UI dialogue
% This opens a figure where you can request setting an arbitrary amount of
% fields (e.g. parameter collection through UI) and then when finished
% collect these into a cell array.
% package: scanpix.helpers
%
% Syntax:  
%   output = scanpix.helpers.makeCustomUIDialogue(prompts)
%   output = scanpix.helpers.makeCustomUIDialogue(prompts, defaultVals)
%
% Inputs: 
%   prompts     - cell array of field name strings 
%   defaultVals - cell array of default values to allocate to prompts (optional)  
%
% Outputs: 
%   output      - cell array (length (propts),2) with field (col1) - value (col2) pairs
%
% See also:
%
% LM 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TO-Do:
% check pix scaling on smaller screen

%% check inputs
if nargin == 1
    defaultVals = cell(length(prompts),1); % all empty fields
end

if length(prompts) ~= length(defaultVals)
    error('Length of field names and default values doesn''t match mate');
end

logInd = cellfun(@islogical,defaultVals);
if any(logInd); defaultVals(logInd) = num2cell(double([defaultVals{logInd}])); end


nLines   = 20;
lineSz1 = [120 20]; % need to check if this scales well on e.g. laptop!
if nargin < 3
    lineSz2 = [120 20];  % need to check if this scales well on e.g. laptop!
else
    lineSz2 = [fillBoxLength 20]; 
end

%% create the dialogue
screenSz = get(0,'screensize');
% dynamically control figure size
if length(prompts) <= nLines
    ySF   = 0.7*(length(prompts)/nLines);
    figSz = [0.425*screenSz(3) (1-ySF)/2*screenSz(4)-25 1.2*(lineSz1(1) + lineSz2(1)) ySF*screenSz(4)+25];
else
    xSF   = 0.2 * ceil(length(prompts)/nLines);
    figSz = [(1-xSF)/2*screenSz(3) 0.175*screenSz(4)-25 xSF*screenSz(3) 0.7*screenSz(4)+25];   % will fail for too many lines, but that seems unlikely? Should maybe include a graceful exit
end

% open figure
fH      = figure('units', 'pixel', 'position', figSz, 'NumberTitle', 'off', 'Name', 'GimmeSomeInput');
figSz   = get(fH,'position');
offSetX = 0.015 * figSz(3);
offSetY = figSz(4) - 1.2*lineSz1(2);

for i = 1:length(prompts)
    
    uicontrol(fH,'units','pixel','position',[ offSetX offSetY lineSz1],'style', 'text', 'string', prompts{i},'horizontalalignment','right','BackgroundColor','W');
    uicontrol(fH,'units','pixel','position',[ offSetX+1.1*lineSz1(1) offSetY+2 lineSz2],'style', 'edit', 'string', defaultVals{i},'horizontalalignment','left');
    
    % update offsets
    offSetY     = offSetY - 0.03*screenSz(4);   
    if offSetY < 0.075*figSz(4)
        offSetX = offSetX + 1.05*offSetX+lineSz1(1)+lineSz2(1);
        offSetY = figSz(4) - 1.2*lineSz1(2);

    end
end
% add button
uicontrol(fH,'Style','PushButton','Units','pixel','Position',[0.4*figSz(3) 0.05*figSz(4) 50 25], 'String', 'Done', 'CallBack', 'uiresume(gcbf)');
% wait for user input
uiwait(fH);

% graceful exit in case dialogue is closed
if ~isgraphics(fH)
   output = {};
   return
end

cnt = 1;
output = prompts(:);
% children are stored bottom to top in parent handle
for i = 2*length(prompts):-2:2
    output{cnt,2} = fH.Children(i).String;
    % now we want to 'guess' numeric input and convert it
    if all(ismember(output{cnt,2}, '0123456789+-.eEdD')) 
        output{cnt,2} = str2num(output{cnt,2}); %#okagrow
    end
    cnt = cnt + 1;
end

if any(logInd); output(logInd,2) = num2cell(logical([output{logInd,2}])); end

close(fH);
end
