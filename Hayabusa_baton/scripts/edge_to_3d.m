function [edge_points_woc, edge_rays] = edge_to_3d(z, fov_angle, trim_u, trim_v)

%Put them in 3D space
% w is going to be the third image dimension

%z = (focal_len*known_width)/app_width;

%info from image's .lbl file
% z = 8.00754;%km
% fov_angle = 5.689; %degrees in vert and horz FOV

horz_dist = z*tand(fov_angle);
vert_dist = z*tand(fov_angle);

%convert pixel numbers to distances in km from origin
dist_u = (trim_u.*(horz_dist/1024));
dist_v = (trim_v.*(vert_dist/1024));

%Make rays
len = length(dist_u(:,1));
dist_w(1:len,1) = -z;


dist_u(len+1,1) = mean(dist_u(:,1));
dist_v(len+1,1) = mean(dist_v(:,1));
dist_w(len+1,1) = mean(dist_w(:,1));

edge_points = [dist_u, dist_v, dist_w];
%nx6 array of end points of line extensions 
%(1:3) are -z direction points, (4:6) are +z direction points
% edge_rays(:,1:3) = edge_points(:,1:3) - [0,0,.5];
% edge_rays(:,4:6) = edge_points(:,1:3) + [0,0,.5];
% 


%%NOTE FOR DEBUGGING: why subtract the last one?
%%do I need the center point? probably yes, that's what makes it
%%bodycentric
%%remember to remove center point after calculating, it could be fucking everything up

for i = 1:len
    edge_points_woc(i,:) = edge_points(i,:) - edge_points(len+1,:);
end

edge_rays(:,1:3) = edge_points_woc(:,1:3) + [0,0,.5];
edge_rays(:,4:6) = edge_points_woc(:,1:3) - [0,0,.5];




end