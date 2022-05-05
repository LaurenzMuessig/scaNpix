function deleteData(obj, type, varargin)
% delete - delete data of individual trials from dacq or npix objects.
% package: scanpix.helpers
%
% Syntax:
%       scanpix.helpers.delete(obj,type)
%       scanpix.helpers.delete(obj,'trials',trialStr)
%       scanpix.helpers.delete(obj,'cells',cellInd)
%
% Inputs:
%    type     - 'trials' or 'cells'
%    varargin:
%             - trialStr - string/cell array of strings; name(s) of trial(s) to be deleted from object
%               if ommited will open UI dialogue to select trial(s)
%             - cellInd - logical or numeric index of cells to be deleted from data
%
% Outputs:
%
% see also:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parse inputs
switch type
    case 'trials'
        if nargin < 3
            [select, loadCheck] = listdlg('PromptString','Select what data to delete:','ListString',obj.trialNames,'ListSize',[160 100]);
            if ~loadCheck
                warning('scaNpix: Deleting data aborted. More data is better anyway.');
                return;
            end
            trialStr = obj.trialNames(select);
        else
            trialStr = varargin{1};
        end
    case 'cells'
        if nargin < 3
%             [select, loadCheck] = listdlg('PromptString','Select what data to delete:','ListString',obj.trialNames,'ListSize',[160 100]);
%             if ~loadCheck
%                 warning('scaNpix: Deleting data aborted. More data is better anyway.');
%                 return;
%             end
%             trialStr = obj.trialNames(select);
        else
            if islogical(varargin{1})
                cellInd = ~varargin{1};
            elseif isnumeric(varargin{1})
                cellInd = true(size(obj.cell_ID,1),1);
                cellInd([varargin{1}]) = false; 
            end
        end
end
    
%% delete
switch type
    case 'trials'
        deleteInd = ismember(obj.trialNames,trialStr);
        fields = fieldnames(obj);
        for i = 1:length(fields)
            if ~any(strcmp(fields{i},obj.fields2spare)) && ~isempty(obj.(fields{i}))
                if isstruct(obj.(fields{i})) && isscalar(obj.(fields{i}))
                    obj.(fields{i}) = structfun(@(x) x(~deleteInd),obj.(fields{i}),'uni',0);
                else
                    obj.(fields{i}) = obj.(fields{i})(~deleteInd);
                end
            end
        end
        
        % remove the last bits if object is now empty
        if isempty(obj.trialNames)
            for i = 1:length(obj.fields2spare)
                if ~strcmp(obj.fields2spare{i},'params')
                    obj.(obj.fields2spare{i}) = [];
                end
            end
            obj.loadFlag  = false;
            scanpix.helpers.changeParams(obj,'default'); % reset defaults
            obj.mapParams = scanpix.maps.defaultParamsRateMaps(class(obj)); % reset defaults
            warning('scaNpix: Back to square one...');
        end
    case 'cells'

        f = fieldnames(obj.spikeData);
        for i = 1:length(f)
            if ~isempty(obj.spikeData.(f{i}){1}) && ~strcmp(f{i},'sampleT')
                obj.spikeData.(f{i}) = cellfun(@(x) x(cellInd),obj.spikeData.(f{i}),'uni',0);
            end
        end
        % need to test whether this works for 2D rate maps as well
        f = fieldnames(obj.maps);
        for i = 1:length(f)
            indEmpty = cellfun('isempty',obj.maps.(f{i}));
            if ~all(indEmpty)
                if all(cellfun(@(x) size(x,1),obj.maps.(f{i})(~indEmpty)) == length(cellInd)) % standard and lin pos maps should be omitted
                    obj.maps.(f{i})(~indEmpty) = cellfun(@(x) x(cellInd,:),obj.maps.(f{i})(~indEmpty),'uni',0);
                end
            end
        end
        obj.cell_ID = obj.cell_ID(cellInd,:);
        if isa(obj,'scanpix.npix')
            obj.cell_Label = obj.cell_Label(cellInd,:);
        end         
end
