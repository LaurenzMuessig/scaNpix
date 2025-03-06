function directions = alignDirections(directions,north,rotDir)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

arguments
    directions {mustBeNumeric}
    north (1,1) {mustBeNumeric} = 270;
    rotDir (1,:)  {mustBeMember(rotDir,{'cw','ccw'})} = 'cw';
end


shiftVal = 360 - north;
directions = mod(directions + shiftVal,360);

if strcmp(rotDir,'ccw')
    directions = directions - 360;
end

end