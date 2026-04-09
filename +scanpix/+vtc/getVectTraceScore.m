function [BFMean,BFSum,TrFMean,TrFSum,BFMask,PoPrFMask,BslFMask] = getVectTraceScore( mapArray,options )
% Version 1 of trace score, find the barrier field simply as a large-enough
% contiguous block of firing.
%
% Input is a cell array of 3 maps, {pre-bsl, probe, post-bsl}

%% Analysis parameters.
arguments
    mapArray {mustBeA(mapArray,'cell')}
    options.bstProbeDef (1,:) {mustBeMember(options.bstProbeDef,{'BFSum','BFMean'})} = 'BFSum';
    options.mapNormVTC (1,:) {mustBeMember(options.mapNormVTC,{'norm2pk','Z'})}  = 'Z';   %   'Z';
    options.BslFieldThr (1,1) {mustBeNumeric} = 0.5;      % 'Pk/2';
    options.BarrFieldThr (1,1) {mustBeNumeric} = 1;      % 'Pk/2';
    options.PostFieldThr (1,1) {mustBeNumeric} = 1;
    options.subBslFiring (1,1) {mustBeNumericOrLogical} = 1;
    options.BFInvalidIfBslFOverlap (1,1) {mustBeNumericOrLogical} = 0;        % If 1, barrF regions which overlap with main field are excluded outright. If 0, main field pixels are still removed from barr F regions, but then these can still count as barr fields.
    options.BFSortParam (1,:) {mustBeMember(options.BFSortParam,{'Area','MeanIntensity'})} = 'Area';   %  'MeanIntensity', 'Area'  % Select 'the' barrier field from many by size or summed rate? (Param names reflect regionprops arguments).
end

%%
[BFMean,BFSum,TrFMean,TrFSum] = deal(nan(size(mapArray,1),1));
[BFMask,PoPrFMask,BslFMask]   = deal( cell(size(mapArray,1),1) );

%%
% normed rate maps
mapNorm = cellfun(@(x) scanpix.vtc.make_zMap( x, options.mapNormVTC ),mapArray,'UniformOutput',false);

%%
for i = 1:size(mapArray,1)
    % If the 'M' input (maps) is empty, it means input is coming from a probe not run. Return the values defined above.
    if isempty( mapArray{i,1} );   continue;    end

       % Define 'baseline' field (main field pre-probe.
    BslFMask{i,1} =  mapNorm{i,1} >=  options.BslFieldThr;

    % If M{2} is empty, this is a cell without a probe - process BslFMask for output, and return (all other outputs are NaN).
    if isempty( mapArray{i,2} )
        BslFMask{i,1}                       = double( BslFMask{i,1} );  % For consisteny, make BslFMask the same format os the other masks,
        BslFMask{i,1}(isnan(mapArray{i,1})) = 3;                   % i.e. double, where 1 is field, 0 no field and 3 unvisited.
        return
    end

    % Subtract pre-bsl from probe and post.
    Ms = { mapNorm{i,2}-mapNorm{i,1},  mapNorm{i,3}-mapNorm{i,1} };

    % Define 'new fields' (new as compared to pre-probe) in cue and post-probe trials.
    bfMaskFin = cell(1,2);
    for itTr=1:2

        % !!!This is the definition of the barrier and post-barrier fields!!!
        %%% TODO - this is assuming firing subtracted, this has to be the case, to remove cells with pre-barr fields, but then what happens to thr?
        if itTr==1
            fMask = Ms{itTr} >= options.BarrFieldThr;
        else
            fMask = Ms{itTr} >= options.PostFieldThr;
        end

        % Remove 'bsl field' bins from barr field mask (as long as we are not just going to exclude them completely, which will happen below).
        if ~options.BFInvalidIfBslFOverlap
            fMask( BslFMask{i,1} ) = false;
        end

        % Label (so as to facilitate bsl F overlap check).
        fLabels = bwlabel( fMask, 4 );

        % If requested, exclude as invalid any contiguous area that overlaps the bsl field.
        if options.BFInvalidIfBslFOverlap
            labelList = setdiff( unique(fLabels), 0 );
            for i=1:length(labelList)
                tempMask = fLabels==labelList(i);
                if any( tempMask(:) & BslFMask{i,1}(:) )
                    fLabels( tempMask ) = 0;
                end
            end
        end

        % If no barrier field, no trace defined. Return an 'all-false' field mask in this case.
        if all( fLabels(:)==0 );  bfMaskFin{itTr}=zeros( size(fMask) );  continue;   end

        %%% TODO: sort by summed firing, not area?
        RP          = regionprops( fLabels, Ms{itTr}, 'Area', 'MeanIntensity' );
        [~,sortInd] = sort( cat(1, RP(1:end).( options.BFSortParam )), 'descend' );
        bfLabelFin  = sortInd(1);

        bfMask                                     = zeros( size(fMask) );  % Default = 0
        bfMask( fLabels==bfLabelFin )              = 1;                     % The actual used mask = 1
        bfMask( fLabels~=bfLabelFin & fLabels~=0 ) = 2;                     % Other fields, but not that chosen as mask are 2.
        bfMask( isnan(mapNorm{i,3}) )              = 3;                     % Non-visited = 3.

        bfMaskFin{itTr}                            = bfMask;
    end
    %
    BFMask{i,1}                      = bfMaskFin{1};
    PoPrFMask{i,1}                   = bfMaskFin{2};
    BslFMask{i,1}                    = double( BslFMask{i,1} );  % For consisteny, make BslFMask the same format os the other two,
    BslFMask{i,1}(isnan(mapNorm{i,3})) = 3;                   % i.e. double, where 1 is field, 0 no field and 3 unvisited.

     % only calculate trace if there is a barrier field.
    if any( BFMask{i,1}(:) )     
        % Barrier field and trace scores are diff in mean Z between pre and post, within BF(probe) mask. Can be mean or summed firing rate.
        BFMean(i,1)  = mean( mapNorm{i,2}( BFMask{i,1}==1 ) );
        BFSum(i,1)   = sum( mapNorm{i,2}(  BFMask{i,1}==1 ) );
        TrFMean(i,1) = mean( mapNorm{i,3}(  BFMask{i,1}==1 ) );
        TrFSum(i,1)  = sum( mapNorm{i,3}(  BFMask{i,1}==1 ) );
        if options.subBslFiring
            BFMean(i,1)  = BFMean(i,1)  - mean( mapNorm{i,1}( BFMask{i,1}==1 ) );
            BFSum(i,1)   = BFSum(i,1)   - sum( mapNorm{i,1}( BFMask{i,1}==1 ) );
            TrFMean(i,1) = TrFMean(i,1) - mean( mapNorm{i,1}(  BFMask{i,1}==1 ) );
            TrFSum(i,1)  = TrFSum(i,1)  - sum( mapNorm{i,1}(  BFMask{i,1}==1 ) );
        end
    end
end

end