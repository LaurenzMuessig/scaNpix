function xyRot = rotatePoints(xy,refPoints)
%UNTITLED3 Summary of this function goes here

% infer from data
% if nargin == 1
%     tmp = xy;
%     tmp(xy(:,2) > 5,:) = NaN;
%     [~,ind] = min(tmp(:,1),[],'omitnan');
%     refPoints = xy(ind,:)';
%     tmp = xy;
%     tmp(xy(:,1) < max(xy(:,1)) - 5,:) = NaN;
%     [~,ind] = min(tmp(:,2));
%     refPoints(:,2) = xy(ind,:)';
% end


% angle of misalignment
theta = atan2(diff(refPoints(2,:)),diff(refPoints(1,:))); 
% rotate
R = [cos(theta) -sin(theta); sin(theta) cos(theta)]; % rotation matrix
xyRot = xy * R;

end

