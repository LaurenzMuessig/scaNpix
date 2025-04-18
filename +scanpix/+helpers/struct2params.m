function prms = struct2params(prmsStruct)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

arguments
    prmsStruct {isstruct(prmsStruct)}
end

prms      = fieldnames(prmsStruct)';
prms(2,:) = struct2cell(prmsStruct)';

end