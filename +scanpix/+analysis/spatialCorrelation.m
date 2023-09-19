function r = spatialCorrelation(mapsA,mapsB)
% spatialCorrelation - Fully vectorised version of computing correlations between rate maps. 
% package: scanpix.analysis
%
% Get the pearson's r by hand so as to vectorise, this is the equation:
% A = A - mean2(A);     B = B - mean2(B);
% r = sum(sum(  A.*B  )) / (  sqrt( sum(sum(  A.*A  )) * sqrt( sum(sum(   B.*B   )) ) );
%
% We'll run all possible correlations.
%
% Note for input: 
% If you want to correlate e.g. cell_1 trialA vs B, cell_2 trial A vs B etc.. 
% then simply format input so trial_A maps = maps(:,1) and those for 
% trial_B = maps(:,2). This will only output the relevant r's instead of
% the full list
%
%
%  Usage:  
%    r = scanpix.analysis.spatialCorrelation(maps)
%
%  Inputs: 
%    mapsA - cell array of rate maps - either {maps,1} or {mapsA,mapsB} 
%           (make sure all maps are of same size)
%
%  Output: 
%    r    - 1:nCells x 1:nCells array of correlations or 1:nCells array of
%           pairwise correlations of maps A v B (see above for input format)
%
% LM 2020


%% some input checks
if ~iscell(mapsA) || (nargin == 2 && ~iscell(mapsB))
   error('Gimme some cell array(s) of maps as input please, will ya?'); 
end

if isempty(mapsA) || (nargin == 2  && isempty(mapsB))
   error('How would I get the correlations for an empty array of maps? Doesn''t work or does it?'); 
end

if nargin == 2
   
    sizeA = length(mapsA);
    sizeB = length(mapsA);
    if ~any(sizeA==1 | sizeB==1) && sizeA ~= sizeB
        error('Dimensions of map arrays are wrong. You can either have a 1x1 map array vs. nx1 map array or an nx1 map array vs. nx1 map array! Now go back and try harder.')
    end
    
    % interpolate maps to common size in case dimensions differ - this is
    % inspired from CBarry's code
    if size(unique([size(mapsA{1});size(mapsB{1})],'rows'),1) > 1
        maxDimSz = max([size(mapsA{1});size(mapsB{1})]);
        for i=1:length(mapsA)
            mapsA{i}=interp2(mapsA{i},linspace(1,size(mapsA{i},2), maxDimSz(2)),linspace(1, size(mapsA{i},1), maxDimSz(1))', '*linear');
        end
        %
        for i=1:length(mapsB)
            mapsB{i}=interp2(mapsB{i},linspace(1,size(mapsB{i},2), maxDimSz(2)),linspace(1, size(mapsB{i},1), maxDimSz(1))', '*linear');
        end
    end
end

%% get the correlations

if nargin == 2
    % Here we run A v B pairwise correlations
    unVisInd = isnan(mapsA{1}) | isnan(mapsB{1});
    unVisIndA = repmat(unVisInd, 1,1, length(mapsA)); % get into correct format;
    unVisIndB = repmat(unVisInd, 1,1, length(mapsB)); % get into correct format;

    rMapArrayA            = cat(3,mapsA{:});
    rMapArrayA(unVisIndA) = NaN;
    A                     = bsxfun(@minus, rMapArrayA , nanmean( nanmean(rMapArrayA ,1), 2 ) );
    %
    rMapArrayB            = cat(3,mapsB{:});
    rMapArrayB(unVisIndB) = NaN;
    B                     = bsxfun(@minus, rMapArrayB , nanmean( nanmean(rMapArrayB ,1), 2 ) );
    %
    AB = bsxfun(@times, A, B);
    sumAB =  nansum(nansum(AB,1),2);                          % This is a 1x1xnCorr vector, one sum(AB) for each corr.
    sqrAA = sqrt( nansum(nansum(A.^2,1),2) );        % This is a 1x1 integer
    sqrBB = sqrt( nansum(nansum(B.^2,1),2) );        % This is a 1x1xnCorr vector, one sqrt(B.^2) for each corr.
    r     = squeeze(sumAB ./ bsxfun(@times, sqrAA, sqrBB));
    
else
    % Here we run all possible correlations in the map array
    fprintf('Running all possible correlation in map array.\n');
    
    rMapArrayA = reshape(horzcat(mapsA{:}), numel(mapsA{1}), length(mapsA));
    rMapArrayA = bsxfun(@minus, rMapArrayA, nanmean(rMapArrayA, 1) ); % subtract mean/cell
    %    
    A         = rMapArrayA(:);                                                     % this is column vector of all rate maps ([rateMap_cell1; rateMap_cell2;...;rateMap_cellN])
    B         = repmat( rMapArrayA, length(mapsA), 1 );                            % these are all ratemap column vectors copied nCell times (row wise), so in the end is ratemap_colVector*nCells x nCells
    AB        = reshape( bsxfun(@times, A, B), numel(mapsA{1}), length(mapsA)^2 ); % need to reshape so each column corresponds to a ratemap col vector for doing the sums
    sumAB     = nansum(AB, 1);                                                     % This is a 1 x nCells^2 vector, one sum(AB) for each ratemap comparison.
    sqrAA     = sqrt( nansum(rMapArrayA.^2, 1) )';                                 % only need to calculate once as product sqrt(AA) x sqrt(BB) == sqrt(AA) x transpose(sqrt(AA))
    sqrAABB   = sqrAA * sqrAA';
    r         = sumAB ./ sqrAABB(:)';                                              % this is a length(maps)^2 row vector of all correlations
    r         = reshape(r, length(mapsA), length(mapsA));                          % reshape into confusion matrix, so rows = 1:cellN vs column 1:cellN ( r(i,j) == r(j,i) )
 
end

end

