function [maps, trimInd] = trimMaps(maps,trimVal,trimIndOnly)

arguments
    maps
    trimVal (1,1) = NaN;
    trimIndOnly (1,1) {mustBeNumericOrLogical} = false;
end
%%
if ~iscell(maps)
    maps = {maps};
    reformat = true;
else
    reformat = false;
end

%%
if isnan(trimVal)
    trimInd(1,1) = find(sum(isnan(maps{1}),2) ~= size(maps{1},2),1,'first');
    trimInd(1,2) = find(sum(isnan(maps{1}),2) ~= size(maps{1},2),1,'last');
    trimInd(2,1) = find(sum(isnan(maps{1}),1) ~= size(maps{1},1),1,'first');
    trimInd(2,2) = find(sum(isnan(maps{1}),1) ~= size(maps{1},1),1,'last');
else
    trimInd(1,1) = find(sum(maps{1} == trimVal,2)~= size(maps{1},2),1,'first');
    trimInd(1,2) = find(sum(maps{1} == trimVal,2)~= size(maps{1},2),1,'last');
    trimInd(2,1) = find(sum(maps{1} == trimVal,1) ~= size(maps{1},1),1,'first');
    trimInd(2,2) = find(sum(maps{1} == trimVal,1) ~= size(maps{1},1),1,'last');
   
end
% trim
if ~trimIndOnly
    maps = cellfun(@(x) x(trimInd(1,1):trimInd(1,2),trimInd(2,1):trimInd(2,2)),maps,'uni',0);
    if reformat; maps = cell2mat(maps); end
end

end