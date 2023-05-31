function r = spatialCorrelation(maps)
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
%    maps - cell array of rate maps - either {maps,1} or {mapsA,mapsB} 
%           (make sure all maps are of same size)
%
%  Output: 
%    r    - 1:nCells x 1:nCells array of correlations or 1:nCells array of
%           correlations of trial A v B (see above for input format)
%
% LM 2020


%% some input checks
if ~iscell(maps)
   error('Gimme some cell array of maps as input please, will ya?'); 
end

if isempty(maps)
   error('How would I get the correlations for an empty array of maps? Doesn''t work or does it?'); 
end


% switch between modes
if size(maps,2) > 1
    % A v B
    maps   = vertcat(maps(:));
    ABflag = true;
    fprintf('Running ''scanpix.analysis.spatialCorrelation'' in A v B mode.\n');
else
    % the full monty
    ABflag = false;
end

% interpolate maps to common size in case dimensions differ - this is
% inspired from CBarry's code
if size( unique( cell2mat( cellfun(@(x) size(x), maps, 'uni', 0) ), 'rows' ), 1) > 1 %
    maxDimSz = max(unique( cell2mat( cellfun(@(x) size(x), maps, 'uni', 0) ), 'rows'));
    for i=1:length(maps)
        maps{i}=interp2(maps{i},linspace(1,size(maps{i},2), maxDimSz(2)),linspace(1, size(maps{i},1), maxDimSz(1))', '*linear');
    end
end

%% get the correlations
% this only does anything if you want to get A v B
unVisInd = isnan(maps{1}) | isnan(maps{round(length(maps)/2)+1});
unVisInd = repmat(unVisInd(:), 1, length(maps)); % get into correct format;

% first arrange ratemaps in a ratemap_colVector x nCells array (so each column is one rate map)
rMapArray           = reshape(horzcat(maps{:}), numel(maps{1}), length(maps));
rMapArray(unVisInd) = NaN; % set unvisited to NaN (only relevant if you compare different trials)
rMapArray           = bsxfun(@minus, rMapArray, nanmean(rMapArray, 1) ); % subtract mean/cell

A         = rMapArray(:);                                                    % this is column vector of all rate maps ([rateMap_cell1; rateMap_cell2;...;rateMap_cellN])
B         = repmat( rMapArray, length(maps), 1 );                            % these are all ratemap column vectors copied nCell times (row wise), so in the end is ratemap_colVector*nCells x nCells
AB        = reshape( bsxfun(@times, A, B), numel(maps{1}), length(maps)^2 ); % need to reshape so each column corresponds to a ratemap col vector for doing the sums
sumAB     = nansum(AB, 1);                                                   % This is a 1 x nCells^2 vector, one sum(AB) for each ratemap comparison.
sqrAA     = sqrt( nansum(rMapArray.^2, 1) )';                                % only need to calculate once as product sqrt(AA) x sqrt(BB) == sqrt(AA) x transpose(sqrt(AA))
sqrAABB   = sqrAA * sqrAA';
r         = sumAB ./ sqrAABB(:)';                                            % this is a length(maps)^2 row vector of all correlations
r         = reshape(r, length(maps), length(maps));                          % reshape into confusion matrix, so rows = 1:cellN vs column 1:cellN ( r(i,j) == r(j,i) )


% A v B between cells (diagonal of upper right quadrant)  
if ABflag
    sz    = length(r)/2;
    r     = diag(r(1:sz,sz+1:end));
end


end

