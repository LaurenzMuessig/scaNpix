function Res = classifySpeedTuning(obj,trialInd,cellInd,options)
%UNTITLED5 Summary of this function goes here

arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    trialInd (1,1) {mustBeNumeric} = 1;
    cellInd {mustBeNumericOrLogical} = true(length(obj.cell_ID(:,1)),1);
    options.minSpeed (1,1) {mustBeNumeric} = 2;
    options.opts (1,1) {mustBeA(options.opts,'optim.options.Fmincon')} = optimoptions('fmincon','Display','off','Algorithm','interior-point');
    options.plot (1,1) {mustBeNumericOrLogical} = false;
    options.add2obj (1,1) {mustBeNumericOrLogical} = false;
end

%%
nParamsUnModel  = 1;
nParamsLinModel = 2;
nParamsSatModel = 3;           

%%
spkTimes = obj.spikeData.spk_Times{trialInd}(cellInd);
%
if isempty(obj.maps.speed{trialInd}); obj.addMaps('speed',trialInd); end
speedMaps = cellfun(@(x) x(:,1:2),obj.maps.speed{trialInd}(cellInd),'UniformOutput',false);
%    
if islogical(cellInd)
    outInd = find(cellInd);
else
    outInd = cellInd;
end

%%
speed     = obj.posData.speed{trialInd};
speedFilt = speed >= options.minSpeed;
speed     = speed(speedFilt);

%% preallocate
Res = struct('rateHat',nan(length(outInd),1),'betaHat',nan(length(outInd),2),'pHat',nan(length(outInd),3),'LL',nan(length(outInd),3),'AIC',nan(length(outInd),3),'BIC',nan(length(outInd),3));

for i = 1:length(spkTimes)

    %
    if isempty(spkTimes{i}); continue; end

    % spike counts per speed bin
    instSpikeCount  = accumarray(ceil(spkTimes{i} .* obj.trialMetaData(trialInd).posFs),1,size(obj.posData.speed{trialInd}));
    instSpikeCount  = instSpikeCount(speedFilt);
    
    %% uniform fit
    rate0 = mean(speedMaps{i}(:,2));     % good initial guess

    lb = 1e-8;

    Res.rateHat(i,1) = fmincon( @(r) uniform_nll(r,instSpikeCount,1/obj.trialMetaData(trialInd).posFs),rate0,[],[],[],[],lb,[],[],options.opts);
    
    %% linear fit
    % guess initial parameters
    beta0 = [mean(speedMaps{i}(:,2))  0];
    % lower bound
    lb    = [1e-8 -Inf];
    %
    Res.betaHat(i,:) = fmincon(@(b) linear_nll(b,speed,instSpikeCount,1/obj.trialMetaData(trialInd).posFs),beta0,[],[],[],[],lb,[],[],options.opts);

    %% saturating exponential fit
    % guess initial parameters
    baseline = min(speedMaps{i}(:,2));
    gain     = max(speedMaps{i}(:,2))-baseline;
    target   = baseline + 0.63*gain;
    idx      = find(speedMaps{i}(:,1) >= target, 1, 'first');

    if ~isempty(idx)
        k = 1/max(speedMaps{i}(idx,1),1);   % crude estimate
    else
        k = 0.05;
    end
    p0   = [baseline,gain,k];
    lb   = [1e-8 0 0]; % lower bound
    %
    Res.pHat(i,:) = fmincon(@(p) satexp_nll(p,speed,instSpikeCount,1/obj.trialMetaData(trialInd).posFs),p0,[],[],[],[],lb,[],[],options.opts);

    %%
    % neg. LL for both fits
    Res.LL(i,1)  = -uniform_nll(Res.rateHat(i,:),instSpikeCount,1/obj.trialMetaData(trialInd).posFs);
    Res.LL(i,2)  = -linear_nll(Res.betaHat(i,:),speed,instSpikeCount,1/obj.trialMetaData(trialInd).posFs);
    Res.LL(i,3)  = -satexp_nll(Res.pHat(i,:),speed,instSpikeCount,1/obj.trialMetaData(trialInd).posFs);
    % Akaike information critirion
    Res.AIC(i,1) = 2*nParamsUnModel  - 2*Res.LL(i,1);
    Res.AIC(i,2) = 2*nParamsLinModel - 2*Res.LL(i,2);
    Res.AIC(i,3) = 2*nParamsSatModel - 2*Res.LL(i,3);
    % Bayesian information critirion
    Res.BIC(i,1) = log(length(instSpikeCount))*nParamsUnModel - 2*Res.LL(i,1);
    Res.BIC(i,2) = log(length(instSpikeCount))*nParamsLinModel - 2*Res.LL(i,2);
    Res.BIC(i,3) = log(length(instSpikeCount))*nParamsSatModel - 2*Res.LL(i,3);

    if options.add2obj
        obj.maps.speed{trialInd}{outInd(i),3}        = nan(6,3);
        obj.maps.speed{trialInd}{outInd(i),3}(1,1)   = Res.rateHat(i);
        obj.maps.speed{trialInd}{outInd(i),3}(2,1:2) = Res.betaHat(i,:);
        obj.maps.speed{trialInd}{outInd(i),3}(3,:)   = Res.pHat(i,:);
        obj.maps.speed{trialInd}{outInd(i),3}(4,:)   = Res.LL(i,:);
        obj.maps.speed{trialInd}{outInd(i),3}(5,:)   = Res.AIC(i,:);
        obj.maps.speed{trialInd}{outInd(i),3}(6,:)   = Res.BIC(i,:);
    end

    %
    if options.plot
        figure;
        % re-bin purely for display
        edges  = options.minSpeed:2:ceil(max(speedMaps{1}(:,1)));
        binIdx = discretize(speed, edges);
        binCtr = edges(1:end-1)+1;
        dispRate = accumarray(binIdx(~isnan(binIdx)), instSpikeCount(~isnan(binIdx)).*obj.trialMetaData(trialInd).posFs, [length(binCtr) 1], @mean, NaN);
        subplot(1,2,1); scatter(binCtr, dispRate, 25, 'k', 'filled');
        %
        [~, minInd] = min(BIC(outInd(i),:));
        hold on;
        vFit = linspace(options.minSpeed, ceil(max(speedMaps{1}(:,1))), 20)';
        if minInd == 2
            plot(vFit, betaHat(outInd(i),1) + betaHat(outInd(i),2)*vFit, 'b-', 'LineWidth', 1.5);
            title(sprintf('Speed tuning: linear'), 'Interpreter', 'none');
        elseif minInd == 3
            plot(vFit, pHat(outInd(i),1)+pHat(outInd(i),2).*(1-exp(-pHat(outInd(i),3).*vFit)), 'c-', 'LineWidth', 1.5);
            title(sprintf('Speed tuning: sat. exponential'), 'Interpreter', 'none');
        else
            plot(vFit,rateHat(outInd(i)) .* ones(length(vFit),1),'k--', 'LineWidth', 1.5);
            title(sprintf('Speed tuning: uniform'), 'Interpreter', 'none');
        end
        set(gca,'xlim',[options.minSpeed,ceil(max(speedMaps{1}(:,1)))]);
        xlabel('Running speed (cm/s)'); ylabel('Firing rate (Hz)');
        axis square
        %
        hold off
        %
        subplot(1,2,2);
        scanpix.plot.plotSpeedMap(obj.maps.speed{trialInd}{outInd(i),1},'ax',gca);
        hold on
        
        if minInd == 2
            plot(betaHat(outInd(i),2).*speedMaps{i}(:,1)+betaHat(outInd(i),1),'b-', 'LineWidth', 1.5);
            title(sprintf('Speed tuning: linear'), 'Interpreter', 'none');
        elseif minInd == 3
            plot(pHat(outInd(i),1)+pHat(outInd(i),2).*(1-exp(-pHat(outInd(i),3).*speedMaps{i}(:,1))),'c', 'LineWidth', 1.5);
            title(sprintf('Speed tuning: sat. exponential'), 'Interpreter', 'none');
        else
            plot(rateHat(outInd(i)) .* ones(length(speedMaps{i}(:,1)),1),'k--', 'LineWidth', 1.5);
            title(sprintf('Speed tuning: uniform'), 'Interpreter', 'none');
        end
        axis square
        hold off
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nll = uniform_nll(rate,spikes,dt)

rate = max(rate,1e-10);

lambda = rate*dt;

nll = -sum(spikes)*log(lambda) + numel(spikes)*lambda + sum(gammaln(spikes+1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nll = linear_nll(beta,speed,spikes,dt)

rate = beta(1) + beta(2)*speed;
% avoid log(0)
rate(rate<1e-10)=1e-10;
lambda          = rate*dt;
%
nll             = -sum(spikes.*log(lambda) - lambda - gammaln(spikes+1)); % gammaln(x+1) == log(n!)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nll = satexp_nll(p,speed,spikes,dt)

rate = p(1) + p(2)*(1-exp(-p(3)*speed));

rate(rate<1e-10)=1e-10;

lambda = rate*dt;

nll = -sum(spikes.*log(lambda)-lambda-gammaln(spikes+1));

end