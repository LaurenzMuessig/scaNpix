function xyRot = rotatePoints(xy,refPoints)
%UNTITLED3 Summary of this function goes here

% angle of misalignment
theta = atan2(diff(refPoints(2,:)),diff(refPoints(1,:))); 
% rotate
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
xyRot = xy * R;

end

