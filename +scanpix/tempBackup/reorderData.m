function reorderData(obj, orderIndex)
% reorderData - reorder trial sequence in dacq or npix class object
% e.g. after adding or removing data
% package: scanpix.helpers
%
% Syntax:
%       scanpix.helpers.reorderData
%       scanpix.helpers.obj.reorderData(orderIndex)
%
% Inputs:
%    orderIndex - numeric index; will order data in object accordingly
%                 (if ommited will open UI dialogue to indicate index)
%
%
% Outputs:
%
% see also:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    uiInput = inputdlg({obj.trialNames{:}},'Please indicate new index for each trial');
    if isempty(uiInput)
        warning('scaNpix: Reordering of trials aborted. Just ask yourself why you started it then...');
        return;
    end
    orderIndex = str2num(cell2mat(uiInput))';
    
    if any(orderIndex > length(obj.trialNames))
        ME = MException('scaNpix:reorderData:orderIndex','OrderIndex for trials > number of trials. That can''t work, can it?');
        throw(ME);
    end
end

% re-order
fields = fieldnames(obj);
for i = 1:length(fields)
    if ~any(strcmp(fields{i},obj.fields2spare)) && ~isempty(obj.(fields{i}))
        if isstruct(obj.(fields{i})) && isscalar(obj.(fields{i}))
            obj.(fields{i}) = structfun(@(x) x(orderIndex),obj.(fields{i}),'uni',0);
        else
            obj.(fields{i}) = obj.(fields{i})(orderIndex);
        end
    end
end

end
