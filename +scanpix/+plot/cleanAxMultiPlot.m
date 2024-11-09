function cleanAxMultiPlot(axArray,mode)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%%
arguments
  axArray cell 
  mode    char {mustBeMember(mode,{'del','inv'})}
end

%%
% empty axes in ax array
remInd = cellfun(@(x) isempty(get(x, 'Children')),axArray);
%
switch mode
    case 'del'
        % delete axes
        delete([axArray{remInd}]);
    case 'inv'
        % make invisible
        cellfun(@(x) axis(x,'off'), axArray(remInd));
end

end