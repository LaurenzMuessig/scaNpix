function output = makeCustomUIDialogue(prompts, defaultVals,fillBoxLength)
% makeCustomUIDialogue - create a custom UI dialogue
% This opens a figure where you can request setting an arbitrary amount of
% fields (e.g. parameter collection through UI) and then when finished
% collect these into a cell array.
% Numeric/logical defaults can be scalars, vectors or matrices; they are
% displayed in MATLAB syntax (e.g. '[1 2 3]') and anything that evaluates
% to a numeric array can be typed in (e.g. '1 2 3', '[1,2;3,4]', '0:5:20').
% package: scanpix.helpers
%
% Syntax:  
%   output = scanpix.helpers.makeCustomUIDialogue(prompts)
%   output = scanpix.helpers.makeCustomUIDialogue(prompts, defaultVals)
%   output = scanpix.helpers.makeCustomUIDialogue(prompts, defaultVals, fillBoxLength)
%
% Inputs: 
%   prompts       - cell array of field name strings 
%   defaultVals   - cell array of default values to allocate to prompts (optional)
%                   Output type follows the default: numeric -> numeric,
%                   logical -> logical, non-empty char -> char, empty ->
%                   numeric if the entry looks numeric, char otherwise
%   fillBoxLength - width (pixel) of edit boxes (optional; default: 120)
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
%%
arguments
    prompts {mustBeA(prompts,'cell')}
    defaultVals {mustBeA(defaultVals,'cell')} = cell(length(prompts),1);
    fillBoxLength (1,1) {mustBeNumeric} = 120;
end

%%
if length(prompts) ~= length(defaultVals)
    error('Length of field names and default values doesn''t match mate');
end

%%
nLines  = 20;
lineSz1 = [120 20]; % need to check if this scales well on e.g. laptop!
lineSz2 = [fillBoxLength 20]; 

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
fH      = figure('units', 'pixel', 'position', figSz, 'NumberTitle', 'off', 'Name', 'GimmeSomeInput','ToolBar','none','MenuBar','none');
figSz   = get(fH,'position');
offSetX = 0.015 * figSz(3);
offSetY = figSz(4) - 1.2*lineSz1(2);

hEdit = gobjects(length(prompts),1);
for i = 1:length(prompts)
    
    uicontrol(fH,'units','pixel','position',[ offSetX offSetY lineSz1],'style', 'text', 'string', prompts{i},'horizontalalignment','right','BackgroundColor','W');
    hEdit(i) = uicontrol(fH,'units','pixel','position',[ offSetX+1.1*lineSz1(1) offSetY+2 lineSz2],'style', 'edit', 'string', val2str(defaultVals{i}),'horizontalalignment','left');
    
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

output = prompts(:);
for i = 1:length(prompts)
    output{i,2} = str2val(strtrim(hEdit(i).String), defaultVals{i}, prompts{i});
end

close(fH);
end

%%
function str = val2str(val)
% convert default value to string for edit box
if isempty(val) && ~ischar(val) && ~isstring(val)
    str = '';
elseif islogical(val)
    str = mat2str(double(val)); % show as 0/1
elseif isnumeric(val)
    str = mat2str(val, 6);
elseif ischar(val) || isstring(val)
    str = char(val);
else
    error('scaNpix::helpers::makeCustomUIDialogue:Can''t display default value of class ''%s''.',class(val));
end
end

%%
function val = str2val(str, defVal, prompt)
% convert edit box string back to value, using default to decide on type
if (isnumeric(defVal) || islogical(defVal)) && ~isempty(defVal)
    if isempty(str)
        val = defVal([]); % keep class
        return
    end
    [val, ok] = parseNumeric(str);
    if ~ok
        warning('scaNpix::helpers::makeCustomUIDialogue:Couldn''t parse ''%s'' for ''%s'' as numeric. Using default value instead.',str,prompt);
        val = defVal;
    elseif islogical(defVal)
        val = logical(val);
    end
elseif ~isempty(defVal)
    % non-empty char default - keep as char
    val = str;
else
    % no default - guess if numeric
    [val, ok] = parseNumeric(str);
    if ~ok; val = str; end
end
end

%%
function [val, ok] = parseNumeric(str)
% parse numeric arrays like '1 2 3', '[1,2;3,4]', '0:5:20', '-Inf NaN 1e3'.
% Only evaluate strings made from safe characters (str2num uses eval)
val = [];
ok  = false;
if isempty(str); return; end
tmp = regexprep(str,'\<(NaN|nan|Inf|inf|pi|true|false)\>','');
if isempty(regexp(tmp,'^[\d\s\.,;:\[\]\+\-\*/\(\)eEdD]*$','once')); return; end
[val, ok] = str2num(str); %#ok<ST2NM>
ok = ok && (isnumeric(val) || islogical(val));
end
