function sparsity = getSparsity(rMap,posMap)
% getSparsity - compute sparsity of rate map (Skaggs et al., 1996)
% package: scanpix.analysis
%
%
% Syntax:
%       sparsity = scanpix.analysis.getSparsity(rMap,posMap)
%
% Inputs:
%    rMap     - rate map
%    posMap   - position map
%
% Outputs: 
%
%    sparsity - sparsity of rate map
%
% LM 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
arguments
    rMap {mustBeNumeric}
    posMap {mustBeNumeric}
end

%%
P_x = posMap ./ sum(posMap(:),'omitnan');
%
num   = sum(rMap(:) .* P_x(:), 'omitnan')^2;
denom = sum(rMap(:).^2 .* P_x(:), 'omitnan');
%
sparsity = num / denom;

end