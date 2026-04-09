function prms = getParamsVTC
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

%%
prms.general.traceFieldOnly = 1;
prms.general.mapNormMode    = 'Z'; 

%%
prms.VTCscore.bstProbeDef            = 'BFSum';
prms.VTCscore.mapNormVTC             = prms.general.mapNormMode ;   %   'Z';
prms.VTCscore.BslFieldThr            = 0.5;      % 'Pk/2';
prms.VTCscore.BarrFieldThr           = 1;      % 'Pk/2';
prms.VTCscore.PostFieldThr           = 1;
prms.VTCscore.subBslFiring           = 1;
prms.VTCscore.BFInvalidIfBslFOverlap = 0;        % If 1, barrF regions which overlap with main field are excluded outright. If 0, main field pixels are still removed from barr F regions, but then these can still count as barr fields.
prms.VTCscore.BFSortParam            = 'Area';   %  'MeanIntensity', 'Area'  % Select 'the' barrier field from many by size or summed rate? (Param names reflect regionprops arguments). 

%%
prms.Zmaps.mode = prms.general.mapNormMode;

%%
prms.VectMaps.spatBinForVMap  = 10; 
prms.VectMaps.wallSegOverSamp = 10;   % Relative to spatial map bins.
prms.VectMaps.vMapBinAng      = 6;    % Units = deg
prms.VectMaps.vMapBinDist     = 1;    % Units = *spatial* map bins. IMPORTANT: see note May 2020 in bvcTrVectMap regarding this parameter (basically must always be 1).
prms.VectMaps.maxEnvSz        = 75;  % Largest box dim = 150cm, limit max tuning to half this. ppm=400.
prms.VectMaps.nBinsEdgeThr    = 15;
prms.VectMaps.excOuterWall    = 1;
prms.VectMaps.smoothVMap      = 1;
prms.VectMaps.smoothKDist     = 5;   % Kernel size for distance, units=vector map bins.
prms.VectMaps.smoothKAng      = 30; 

%%
prms.plotVectMaps.nSteps    = 10;
prms.plotVectMaps.isCirc    = true;
prms.plotVectMaps.vMapsNCrc = 4;

%%
prms.fMasks.maskCM = [0.847058824 0.847058824 0.847058824; 0 0.439215686 0.749019608; 1 0.937254902 0; 0.568627451 0.819607843 0.309803922;1 1 1];

end