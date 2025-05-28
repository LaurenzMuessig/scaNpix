function [dataIndex, missingTrials] = matchTrialSeq2Pattern(trialSeqObj,pattern2extract,options)  
% matchTrialSeq2Pattern - match a pattern of trial types to the trials run for a specific dataset 
% package: scanpix.helpers
%
% Syntax:  
%    scanpix.helpers.matchTrialSeq2Pattern(trialSeqObj,pattern2extract)
%
% Inputs:
%    trialSeqObj     - cell array of trial type sequence
%    pattern2extract - cell array of trial types to extract 
%    varargin        -  
%
%
% Outputs:
%    dataIndex     - numeric index to extract trials in order of 'pattern2extract'
%    missingTrials - numeric index of missing trial types in 'trialSeqObj'
%
% LM 2024
%
%% 
arguments
  trialSeqObj (1,:) {mustBeA(trialSeqObj,'cell')}
  pattern2extract (1,:) {mustBeA(pattern2extract,'cell')}
  options.mode (1,:)  {mustBeMember(options.mode,{'pattern','exact'})} = 'pattern';
  options.bslKey (1,:) {mustBeA(options.bslKey,'cell')} = {'fam'};
  options.ignKey (1,:) {mustBeA(options.ignKey,'cell')} = {'sleep'};
  options.getFlankBSL (1,1) {mustBeNumericOrLogical} = true;
  options.exactflag (1,1) {mustBeNumericOrLogical} = true;
end

%%
switch options.mode
    case 'pattern'
        % 
        if ~isempty(pattern2extract{2})
            if options.exactflag
                checkValid = all(ismember(unique(pattern2extract{2}),trialSeqObj));
            else
                uniqueTrialsPattern = unique(pattern2extract{2});
                checkValid = true;
                for i = 1:length(uniqueTrialsPattern)
                    checkValid = checkValid & any(~cellfun('isempty',regexp(trialSeqObj,uniqueTrialsPattern{i},'once')));
                end
            end
        else
            checkValid = true;
        end
        %
        if ~checkValid
            dataIndex     = [];
            missingTrials = 1:length([pattern2extract{:}]);
            return;
        end

        %
        dataIndex           = zeros(2,length([pattern2extract{:}])+2*options.getFlankBSL);
        nPreProbeBSLPattern = length(pattern2extract{1});

        firstNonBslInTrialSeq           = find( ~ismember(trialSeqObj,[options.bslKey,options.ignKey]), 1, 'first' );
        if isempty(firstNonBslInTrialSeq)
            firstNonBslInTrialSeq = length(trialSeqObj) + 1;
        end

        blsIndex                        = find(strcmp(trialSeqObj,unique(pattern2extract{1})),nPreProbeBSLPattern,'first');
        % if ~(all(blsIndex < firstNonBslInTrialSeq) || all(blsIndex > firstNonBslInTrialSeq))
        %     blsIndex = blsIndex(blsIndex < firstNonBslInTrialSeq);
        % end
        dataIndex(1,1:length(blsIndex)) = blsIndex;
        dataIndex(2,1:length(blsIndex)) = blsIndex < firstNonBslInTrialSeq;
        
        % graceful exit
        if isempty(pattern2extract{2}) || isempty(blsIndex)
            missingTrials               = find(dataIndex(1,:)==0);
            dataIndex(:,missingTrials)  = [];
            return
        end

        % now do all trials after the (pre-probe) baseline trials
        postBSLTrial2Extr    = pattern2extract{2};
        % first check if all pre-probe baselines are correct (in case there was only 1 pre-probe baseline)
        firstProbeInTrialSeq = find( ismember(trialSeqObj,postBSLTrial2Extr), 1, 'first' );
        dataIndex(:,dataIndex(1,:) > firstProbeInTrialSeq) = 0;
        %
        remainTrials      = trialSeqObj(max(dataIndex(1,:))+1:end);
        %
        c1 = max(dataIndex(1,:));
        c2 = nPreProbeBSLPattern + 1 + options.getFlankBSL;
        while ~isempty([postBSLTrial2Extr{:}])

            % partial matching by common string also possible
            if options.exactflag
                tempInd = find( ismember( remainTrials,postBSLTrial2Extr{1} ), 1, 'first' );
            else
                tempInd = find(~cellfun('isempty',regexp(remainTrials,postBSLTrial2Extr{1},'once')), 1, 'first');
            end
            %
            if isempty(tempInd)
                postBSLTrial2Extr = postBSLTrial2Extr(2:end);
            else
                dataIndex(1,c2)       = tempInd + c1;
                postBSLTrial2Extr     = postBSLTrial2Extr(2:end);
                remainTrials(tempInd) = [];
                c1                    = c1 + 1;
            end
            c2                        = c2+1;
        end
        % add baselines flanking the probe(s) if desired
        if options.getFlankBSL
            preProbeExtraBSL = find(strcmp(trialSeqObj(max(dataIndex(1,1:length(pattern2extract{1})),[],'omitnan')+1:max(dataIndex(1,:),[],'omitnan')),options.bslKey),1,'last');
            if ~isempty(preProbeExtraBSL); dataIndex(1,length(pattern2extract{1})+1) = preProbeExtraBSL + max(dataIndex(1,1:length(pattern2extract{1})),[],'omitnan'); end
            %
            postProbeExtraBSL = find(strcmp(trialSeqObj(max(dataIndex(1,:),[],'omitnan')+1:end),options.bslKey),1,'first');
             if ~isempty(postProbeExtraBSL); dataIndex(1,end) = postProbeExtraBSL + max(dataIndex(1,:),[],'omitnan'); end
        end
        %
        missingTrials              = find(dataIndex(1,:)==0);
        dataIndex(:,missingTrials) = [];
    
    % match sequence exactly
    case 'exact'
        % graceful exit
        if options.exactflag && ~all(ismember( unique(pattern2extract),trialSeqObj ))
            dataIndex     = [];
            missingTrials = 1:length(pattern2extract);
            return;
        end
        
        % this ends up overly complicated but I can't find an easier way to
        % deal with repeating trial types and differences in the sequence
        % of trials
        % dataIndex = zeros(2,length(pattern2extract));
        % minMax    = nan(length(pattern2extract),2);
        % tempInd   = cell(length(pattern2extract),1);

        % first get the indices for all trial types
        for i = 1:length(trialSeqObj)
            if options.exactflag
                tempInd = find(ismember( pattern2extract, trialSeqObj{i}));
            else
                tempInd = find(~cellfun('isempty',regexp(pattern2extract,trialSeqObj{i},'once')));
            end
            % this wouldn't account for a case where e.g. only a post-probe
            % baseline would be there
            if i == 1
                dataIndex = tempInd(1);
            else
                dataIndex = [dataIndex, tempInd(find(tempInd > max(dataIndex),1,'first'))];
            end
        end 
        dataIndex(2,:) = 0;
        %
        missingTrials  = find(~ismember(1:length(pattern2extract),dataIndex(1,:)));
end

end

%   dataIndex = nan(1,length(pattern2extract));
% 
%         find( ~strcmp(pattern2extract,p.Results.bslk), 1, 'first' );
% 
%         firstNonBslInPattern  = find( ~strcmp(pattern2extract,p.Results.bslk), 1, 'first' );
%         if isempty(firstNonBslInPattern)
%             nBSL = length(pattern2extract);
%         else
%             nBSL = firstNonBslInPattern - 1;
%         end
% 
% 
%         %
%         % if p.Results.ignoreOrder
%         %     bslInd = find( strcmp(trialSeqObj,p.Results.bslk));
%         %     nBSL   = firstNonBslInPattern - 1;
%         %     if length(bslInd) > nBSL
%         %         % bslInd = bslInd(end-(nBSL-1):end);
%         %         bslInd = bslInd(1:nBSL);
%         %     end
%         % else
% 
%         % end
% 
%         if isempty(firstNonBslInTrialSeq)
%             dataIndex     = find(strcmp(trialSeqObj,p.Results.bslk),nBSL,'last'); % extract the n last baseline trials
%             if length(dataIndex) < nBSL
%                 missingTrials = length(dataIndex)+1:nBSL;
%             else
%                 missingTrials = [];
%             end
%             return
%         end
% 
% 
% % % now do all trials after the (pre-probe) baseline trials
% postBSLTrial2Extr = pattern2extract{2};
% if p.Results.ignoreOrder
%     remainTrials      = trialSeqObj(~strcmp(trialSeqObj,p.Results.bslk));
% else
%     remainTrials      = trialSeqObj(firstNonBslInTrialSeq:end);
% end
% %
% c1 = firstNonBslInTrialSeq - 1;
% c2 = nBSL + 1;
% while ~isempty([postBSLTrial2Extr{:}])
% 
%     % partial matching by common string also possible
%     if p.Results.exact
%         tempInd = find( strcmp( remainTrials,postBSLTrial2Extr{1} ), 1, 'first' );
%     else
%         tempInd = find(~cellfun('isempty',regexp(remainTrials,postBSLTrial2Extr{1},'once')));
%     end
%     %
%     if isempty(tempInd)
%         postBSLTrial2Extr = postBSLTrial2Extr(2:end);
%     else
%         dataIndex(c2)     = tempInd + c1;
%         postBSLTrial2Extr = postBSLTrial2Extr(2:end);
%         remainTrials      = remainTrials(2:end);
%         c1                = c1 + 1;
%     end
%     c2                    = c2+1;
% end
% %
% missingTrials             = find(isnan(dataIndex));
% dataIndex(missingTrials)  = [];

% dataIndex = nan(1,length(pattern2extract));
% % Get the n last bsl trials before the first probe.
% if p.Results.ignoreOrder
%     dataIndex(1:length(bslInd)) = bslInd;
% else
%     preProbeBSLInd                      =  find( strcmp(trialSeqObj,p.Results.bslk), nBSL, 'first' );
%     dataIndex(1:length(preProbeBSLInd)) =  preProbeBSLInd;
% 
%     % for i = 1:nBSL
%         % if firstNonBslInTrialSeq-i < 1
%         %     continue
%         % else
%         %     dataIndex( nBSL - (i-1) ) = firstNonBslInTrialSeq - i;
%         % end
% 
%     % end
% end

 % % firstBslInTrialSeq    = find( ~ismember(trialSeqObj,[options.bslKey,options.ignKey]), 1, 'first' );
        % if isempty(firstNonBslInTrialSeq)
        %     blsIndex    = find(strcmp(trialSeqObj,unique(pattern2extract{1})),nPreProbeBSLPattern,'first'); % extract the n last baseline trials
        %     remTrialInd = max(blsIndex) + 1;
        % else
        %     blsIndex    = find(strcmp(trialSeqObj(1:firstNonBslInTrialSeq),unique(pattern2extract{1})),nPreProbeBSLPattern,'first'); % extract the n last baseline trials
        %     remTrialInd = firstNonBslInTrialSeq;
        % end
        % dataIndex(1:length(blsIndex)) = blsIndex;



        % for i = 1:length(pattern2extract)
        %     if options.exactflag
        %         tempInd{i} = find(ismember( trialSeqObj,pattern2extract{i}));
        %     else
        %         tempInd{i} = find(~cellfun('isempty',regexp(trialSeqObj,pattern2extract{i},'once')));
        %     end
        % 
        %     %            
        %     minMax(i,1) = min(tempInd{i});
        %     minMax(i,2) = max(tempInd{i} );
        %     if i > 1
        %         if all(minMax(i,2) < max(minMax(1:i-1,:),[],'omitnan')) || all(minMax(i,2) == max(minMax(1:i-1,:),[],'omitnan'))
        %             minMax(i,:) = NaN;
        %         end
        %     end
        % end
        % 
        % remInd           = all(isnan(minMax),2);
        % indBump          = cumsum(remInd);
        % minMax(remInd,:) = [];
        % tempInd(remInd)  = [];
        % indBump(remInd)  = [];
        % % now extract the actual trial indices
        % for i = 1:size(minMax,1)
        % 
        %     if i == 1
        %         if size(minMax,1) == 1
        %             dataIndex(1,i+indBump(i)) = tempInd{i};
        %         else
        %             selInd = find(tempInd{i} <= min(unique(minMax(i+1:end,:)),[],'omitnan'),1,'last');
        %             if ~isempty(selInd)
        %                 dataIndex(1,i+indBump(i)) = tempInd{i}(selInd);
        %             end
        %         end
        %     % 
        %     elseif i == size(minMax,1)
        %         selInd = find(tempInd{i} > dataIndex(1,i-1),1,'first');
        %         if ~isempty(selInd)
        %             dataIndex(1,i+indBump(i)) = tempInd{i}(selInd);
        %         end
        %     else
        %         selInd = find(tempInd{i} > max(dataIndex(1,1:i-1),[],'omitnan') & tempInd{i} <= max(minMax(i+1,:),[],'omitnan'),1,'first');
        %         if ~isempty(selInd) && ~any(tempInd{i}(selInd) > cellfun(@(x) max(x,[],'omitnan'),tempInd(i+1:end)))
        %             dataIndex(1,i+indBump(i)) = tempInd{i}(selInd);
        %         end
        %     end
        % end