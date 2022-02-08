%2/1/2022
%Dahlia Baker
%OREX run file


clear all
close all
clc

load_script = '../../orex_working/orex_test_data.mat';
img_path = ["orex/To_Send"];
ext = 1; %length of ray extension
limb = 1; %1 for include terminator, 0 for limb-only


load(load_script);
addpath('/Users/dahliabaker/Documents/GradSchool/Research/LimbTracing/limbtracingcode/journal_work/limbtracing_debug/scripts/');

for i = 1:length(img_path)
    addpath(img_path(i));   
end
IR = 0;
IRdat = [];
density = 5;
body = 'bennu';
ir_list = [];
%%
[limb_starts, limb_ends, edge_points_bc] = image_to_limbs_orex(names(1:2:end), r(1:2:end), fov_angle,CB(:,:,1:2:end),sun_pos(1:2:end,:),phase(1:2:end),limb,ext,IR,IRdat,density,body,ir_list);


%%
figure
hold on
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')
color = ['r','b','c','m','y','k'];
color = [color, color, color, color, color, color, color, color, color, color, color, color];
for j = 1:35
    %uncomment below to also plot ray lines with the shape points
%      for i = 1:36
%          line([limb_starts{j}(i,1) limb_ends{j}(i,1)], [limb_starts{j}(i,2) limb_ends{j}(i,2)],[limb_starts{j}(i,3) limb_ends{j}(i,3)],'Color', [1,0,0,0.2],'LineWidth',1);
%      end
     scatter3(edge_points_bc{j}(:,1),edge_points_bc{j}(:,2),edge_points_bc{j}(:,3),'filled')
     drawnow;
end
axis('equal')

%%
[shapeEndPnts,shapePnts, shapePntNhats] = shape_from_limbs_orex(limb_ends,limb_starts, r(1:2:end), 10,2,0,[0,0,-1]);

%plots the shapePnt results of shape_from_limbs
ptCloud = pointCloud(shapePnts,'Normal',shapePntNhats);
figure()
pcshow(ptCloud)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')