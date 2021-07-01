%Run Script
%6/29/2021
%Dahlia Baker

%set initial values
load_script = 'mv_bennu/mv_bennu.mat';
img_path = ["mv_bennu/bennu_automated_images_1";"mv_bennu/bennu_automated_images_2"];
fov_angle = 0.8;
phase = 0; %phase of test case - can be an
ext = 1; %length of ray extension
limb = 1; %1 for include terminator, 0 for limb-only


%%
load(load_script);
addpath('scripts/');

for i = 1:length(img_path)
    addpath(img_path(i));   
end
    
[limb_starts, limb_ends, edge_points_bc] = image_to_limbs_orex(img_name, r, fov_angle,CB,sun_pos,phase,limb,ext);

figure
hold on
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')
color = ['r','b','c','m','y','k'];
color = [color, color, color, color, color, color, color, color, color, color, color, color];
for j = 1:length(r)
    %uncomment below to also plot ray lines with the shape points
%      for i = 1:36
%          line([limb_starts{j}(i,1) limb_ends{j}(i,1)], [limb_starts{j}(i,2) limb_ends{j}(i,2)],[limb_starts{j}(i,3) limb_ends{j}(i,3)],'Color', [1,0,0,0.2],'LineWidth',1);
%      end
     scatter3(edge_points_bc{j}(:,1),edge_points_bc{j}(:,2),edge_points_bc{j}(:,3),'filled')
     drawnow;
end
axis('equal')

[shapeEndPnts,shapePnts, shapePntNhats] = shape_from_limbs_orex(limb_starts,limb_ends, r, 10,2,0,[0,0,-1]);

%plots the shapePnt results of shape_from_limbs
ptCloud = pointCloud(shapePnts,'Normal',shapePntNhats);
figure()
pcshow(ptCloud)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')