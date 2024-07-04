function sttc = getSTTC(st1,st2,trialDur,dt)
%UNTITLED5 Summary of this function goes here
%   https://www.jneurosci.org/content/34/43/14288

%%
% default delta
if nargin == 3
    dt = 0.5;
end
%

if isempty(st1) || isempty(st2)
    sttc = NaN;
    return
end


%%
Ta = getT(st1,dt,trialDur);
Tb = getT(st2,dt,trialDur);
%%
P  = getP(st1,st2,dt);
Pa = P / length(st1);
P  = getP(st2,st1,dt);
Pb = P / length(st2);
%%
sttc = 0.5 * ( (Pa-Tb) / (1 - Pa * Tb) + (Pb - Ta) / (1 - Pb * Ta) );

end

function T = getT(st,dt,trialDur)
%

%%
nSpikes = length(st);
maxPossTime = 2 * nSpikes * dt; % max possible time

if nSpikes == 1
    T = maxPossTime;
    % special case of 1 spike
    if st < dt
        T = T - (dt - st); 
    elseif st+dt > trialDur
        T = T - (dt - (trialDur - st));
    end
else
    % remove times of spike overlap
    stDiff    = diff(st);
    stOverlap = stDiff(stDiff < 2*dt);
    T         = maxPossTime - (2 * dt * length(stOverlap) - sum(stOverlap));
    % adjust first and last window if necessary
    if st(1) < dt
        T = T - (dt - st(1));
    end
    if trialDur - st(end) < dt
        T = T - (dt - (trialDur - st(end)));
    end
end
T = T / trialDur;

end

function P = getP(st1,st2,dt)
%

%%
if length(st2) == 1
    ind    = 1;
else
    ind    = interp1(st2, 1:numel(st2), st1, 'nearest', 'extrap');  
end
%
indOverlap = abs(st1 - st2(ind)) < dt;
P          = nnz(indOverlap);


% indLeft = ind - 1;
% indLeft(indLeft < 1) = 1;
% left_DT = abs(st2(indLeft) - st1) < dt;
% %
% right_DT = abs(st2(ind) - st1) < dt;
% inDT = left_DT | right_DT;
% P = nnz(inDT);



end

