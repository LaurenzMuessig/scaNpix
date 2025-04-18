function spikeProps = getWaveFormProps(obj, sampleRate, options)
% getWaveFormProps - Characterise waveforms 
% ATM only works for dacq class objects
% Calculate mean waveform/channel for each cell as well as a few properties (spike width (peak-trough), mean rate and mean 
% of autocorrelation. We are assuming waveforms are as per DACQ native format 'flipped', so peak is positive!
% This is a short version of 'spk_characterisewaveform' by TW
% Note: It is somewhat time consuming to make all autocorrelograms for a large dataset on the fly, so you can store the properties 
% as ['spikeProps_' obj.trialName '.mat'] to disk (setting prms.saveFlag=true). This would then be loaded next time you load the data 
% into an object and run obj.getWaveFormProps (this makes particular sense if your machine is slow or if you are running a big analysis pipeline
%
% Syntax:  spikeProps = getWaveFormProps(obj)
%          spikeProps = getWaveFormProps(obj, optionalInputStruct );
%          spikeProps = getWaveFormProps(obj, 'inputName', inputVal, .. etc .. );
%
% Inputs:
%    varargin - optional prmsStruct or Name-Value pair list to change default values for params 
%
% Options:
%
%   'Fs',          4800  - sample rate for spike channels (default=48kHz)
%   'acWin',       50    - window size for autocorrelograms in ms (default=50ms)
%   'acBin',       1     - bin size autocorrelograms in ms (default=1ms)
%   'saveFlag',    0     - save flag to write data to disk (default=false)
%   'loadFlag',    1     - load flag to load data from disk (default=true) 
%
% Outputs:
%
%   spikeProps
%
%
% See also: ;
%
% LM 2020

%% TO-DO:
% a) haven't added a lot of the stuff from TW's original fucntion as
% never used, but width at half height might be useful addition?
% b) minimum spike n for prop calculation

% prms.load = 0;     % load flag to load data from disk - making a lot of 
                       % ACs is time consuming, so for big analysis runs it 
                       % makes sense to write these property files to disk first
%%
arguments
    obj {mustBeA(obj,'scanpix.ephys')}
    sampleRate (1,1) {mustBeNumeric}
    options.acWin (1,1) {mustBeNumeric} = 50;
    options.acBin (1,1) {mustBeNumeric} = 1; % in degrees
    options.save (1,1) {mustBeNumericOrLogical} = false; 
end

%% DO THAT THANG

[maxAmp,spkWidth,firstMomentAC,meanRate] = deal(cell(length(obj.trialNames),1));
spikeProps                               = nan(length(obj.cell_ID), 6, length(obj.trialNames));
for i = 1:length(obj.trialNames)
    
    % pre-allocate
    firstMomentAC{i}   = nan(length(obj.cell_ID), 1);
    meanRate{i}        = nan(length(obj.cell_ID), 1);
    maxAmp{i}          = nan(length(obj.cell_ID), 1);
    spkWidth{i}        = nan(length(obj.cell_ID), 1);

    % check if there are spike props files so we can skip calculating all
    % properties
    % if prms.loadFlag
    %     try
    %         tmp = load([obj.dataPath 'spikeProps_' obj.trialNames{i} '.mat']);
    %         % we need to protect here against case where we have removed
    %         % cells from object but spikeprops contains full set
    %         f = fieldnames(tmp);
    %         if ~all( ismember(tmp.f{1}(:,1:2), obj.cell_ID(:,1:2), 'rows') )
    %             cellInd = ismember(tmp.f{1}(:,1:2), obj.cell_ID(:,1:2), 'rows');
    %         else
    %             cellInd = true(size(obj.cell_ID,1), 1);
    %         end
    %         % assign data
    %         spikeProps(:,:,i) = tmp.f{1}(cellInd,:);
    %         continue
    %     catch
    %         warning(['scaNpix::analysis::getWaveFormProps: Can''t find spikeProps file for ' obj.trialNames{i} ' in ' obj.dataPath{i} '.']);
    %     end
    % end
    % loop though cells in data set
    for j = 1:length(obj.cell_ID)
        
                
        meanWFs  = squeeze( mean( obj.spikeData.spk_waveforms{i}{j}, 1, 'omitnan') ); % mean waveforms
        
        % Get max2min amplitude, and the maximum amplitude channel %
        [~, maxCh]                       = max( max(meanWFs, [], 1) - min(meanWFs, [], 1) );
        maxAmpWF                         = meanWFs(:,maxCh);
        % find peak
        localMaxInd                      = find( islocalmax(maxAmpWF) );
        [peakAmp,tempInd]                = max( maxAmpWF(localMaxInd) );
        peakInd                          = localMaxInd(tempInd);
        % find trough
        localMinInd                      = find( islocalmin(maxAmpWF(peakInd:end)) ); % find only after peak
        [troughAmp,tempInd]              = min( maxAmpWF(localMinInd + peakInd-1) );
        troughInd                        = localMinInd(tempInd) + peakInd - 1;
        
        % For very bad waveforms, can be no local max or mins on maxAmp channel. In this case, give up and return all as NaN. %
        if isempty(peakAmp) || isempty(troughAmp)
            [maxAmp{i}(j),spkWidth{i}(j),firstMomentAC{i}(j),meanRate{i}(j) ] = deal(NaN);
            warning(['scaNpix::analysis::getWaveFormProps: Waveform for cell' num2str(obj.cell_ID(j,1)) '_tet' num2str(obj.cell_ID(j,2)) ' in ' obj.trialNames{i} ' is bad ' ...
                       '(nSpikes = ' num2str( length(obj.spikeData.spk_Times{i}{j})  ) '). No local max or min found. You might want to remove that shit!']);
            continue;
        else
            maxAmp{i}(j)    = peakAmp - troughAmp;
            spkWidth{i}(j)  = (troughInd - peakInd) * (1/sampleRate*1000); % in ms
        end
         
        % additional params
        firstMomentAC{i}(j) = scanpix.analysis.get1stMomentAC(obj.spikeData.spk_Times{i}{j},options.acWin,options.acBin,obj.trialMetaData(i).duration);
        meanRate{i}(j)      = length(obj.spikeData.spk_Times{i}{j}) / obj.trialMetaData(i).duration;
        
    end
    
    spikeProps(:,1:2,i) = obj.cell_ID(:,1:2);
    spikeProps(:,3,i)   = spkWidth{i};
    spikeProps(:,4,i)   = firstMomentAC{i};
    spikeProps(:,5,i)   = meanRate{i};
    spikeProps(:,6,i)   = maxAmp{i};
    
    % save spikeprops file
    if options.saveFlag
        tmp = squeeze(spikeProps(:,:,i));
        save([obj.dataPath{i} 'spikeProps_' obj.trialNames{i}],'tmp');
    end
    
end

end




