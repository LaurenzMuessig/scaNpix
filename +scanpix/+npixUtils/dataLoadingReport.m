function dataLoadingReport(nFramesBonsai,nFramesSync,BonsaiCorruptFlag)
% 


% this needs a bit more experimentation if we accounted for all eventualities
% in an acceptable manner, esp. for corrupt data
if ~BonsaiCorruptFlag && nFramesBonsai ~= nFramesSync
    if nFramesBonsai - nFramesSync == -1
        % this case isn't 100% clear as frame could be missing anywhere - don't think this ever happens anymore since using CatGT to extract sync pulse times
        disp('Warning. Missmatch between n of pos samples and n of TTLs. -1 frame in pos data, so we assume the last frame in neuropix stream is incomplete and will be removed from the data.');
    elseif nFramesBonsai - nFramesSync == +1
        % this case should be clear and essentially no frame is missing!
        disp('Warning. Missmatch between n of pos samples and n of TTLs. +1 frame in pos data, so the last TTL was deemed incomplete. We will remove last pos frame.');
    elseif nFramesBonsai - nFramesSync > 1
        % this case should happen when animal unplugs during recording
        fprintf('Warning. There are %i more frames in tracking stream compared to neuropix data - assuming that the animal unplugged in recording. If not, you are in trouble\n', nFramesBonsai-nFramesSync);
    else
        fprintf('Warning. Missmatch between n of pos samples and n of TTLs in neuropix data - %i vs. %i. Better go and check out why.\n', nFramesBonsai, nFramesSync);
    end
elseif BonsaiCorruptFlag
    % Note: It took me a while to figure out that FlyCap is sensitive to how
    % it's closed when switching between cameras or changing settings. This
    % results in the camera not sending it's metadata over (i.e. frame
    % count and time stamp). When we haven't got the framecount we can't
    % reconstruct which frames are missing. So by assuming we have an even
    % sampling interval we will introduce some jitter, so at some point the
    % 2 streams will be more and more out of sync (also depending on when 
    % the frames were skipped, the later the better). So for these cases we
    % should check the data, i.e. e.g. do we see an obvious drop in the 
    % spatial properties. 
    if nFramesBonsai == nFramesSync
        % best case scenario - no frames are missing
        disp('Corrupt Point Grey MetaData Logging: No frame N mismatch with neuropix stream! All good. Phew...') 
    elseif nFramesBonsai - nFramesSync == -1
        disp('Corrupt Point Grey MetaData Logging: -1 frame in position data. This shouldn''t really happen my friend. Time for a deep dive...');
    elseif nFramesBonsai - nFramesSync == +1
        disp('Corrupt Point Grey MetaData Logging: +1 frame in position data. Assuming last frame is missing in neuropix data. This should be equal to all is good and you can stop stressing')
    else
        fprintf('Corrupt Point Grey MetaData Logging: Frame N mismatch is %i (Bonsai) vs. %i (Neuropix). We''ll assume even sampling between frames, but you should check carefully if n of mismatch is too big!\n', nFramesBonsai, nFramesSync);
    end
end

end