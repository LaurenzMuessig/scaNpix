function timesSec = convertPointGreyCamTimeStamps(timeStamps)
% converts metadata timestamps recorded in bonsai from point grey camera into seconds
% package: scanpix.npixUtils

% input:     timeStamps: metadata timestamps recorded in bonsai
% output:    times: seconds, normalized such that first time is 0
%
% written by Richard W Allen
% https://groups.google.com/g/bonsai-users/c/WD6mV94KAQs

% extract first cycle (first 7 bits)
cycle1 = bitshift(timeStamps, -25);

% extract second cycle (following 13 bits)
cycle2 = bitand(bitshift(timeStamps, -12), 8191) / 8000; % 8191 masks all but but last 13 bits

% account for overflows in counter
timesSec = cycle1 + cycle2;
overflows = [0; diff(timesSec) < 0];
timesSec = timesSec + (cumsum(overflows) * 128);

% offset such that first time is zero
timesSec = timesSec - min(timesSec);

end
