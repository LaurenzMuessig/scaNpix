function ind = getIndicesFromStartEnd(startInd,endInd)
% Generate a full list of numeric indices from a list of tart and end indices

%%
arguments
    startInd (:,1) {mustBeNumeric}
    endInd (:,1) {mustBeNumeric}
end

if length(startInd) ~= length(endInd); error('scaNpix::helpers::getIndicesFromStartEnd: You need to have equal length lists of start and end indices! Pretty obvious if you ask me...'); end
%%
L               = length(startInd);
F               = cumsum(endInd-startInd+1);
ind             = ones(1,F(end));
ind(1)          = startInd(1);
ind(1+F(1:L-1)) = startInd(2:L)-endInd(1:L-1);
ind             = cumsum(ind);

end