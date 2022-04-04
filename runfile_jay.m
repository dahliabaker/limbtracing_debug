%Run Script
%6/29/2021
%Dahlia Baker
clear all
close all
clc
tic
%set initial values
datasetPath = "new_datasets/90bennu_sim/";img_path = datasetPath+"/bennu_automated_images";body = "bennu";
% datasetPath = "new_datasets/45itokawa_sim/";img_path = datasetPath+"/itokawa_automated_images";body = "itokawa";
fov_angle = 0.8;
ext = 1; %length of ray extension
limb = 0; %1 for include terminator, 0 for limb-only
IR = 0;%0 if no IR data, 1 if yes
density = 5; %number of degrees between silhouette sampled points

if IR == 0
    ir_imgs = [];
else
    if contains(datasetPath,'bennu')
        load('ir_list_b.mat')
        body = body+"_IR";
    elseif contains(datasetPath,'itokawa')
        load('ir_list_i.mat')
        body = body+"_IR";
    end
end
    
%% run image to limbs
load(datasetPath+'bennu_72.mat');
addpath('scripts_limbonly/');

for i = 1:length(img_path)
    addpath(img_path(i));   
end

% [limb_starts, limb_ends, edge_points_bc] = image_to_limbs_lo(img_name(good), r(good), fov_angle,CB(:,:,good),sun_pos(good,:),phase,ext);
[limb_starts, limb_ends, edge_points_bc] = image_to_limbs_lo(img_name,ir_imgs, r, fov_angle,CB,sun_pos,ext);
%% plot raw pixels and rays
figure(4),clf
hold on
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')
color = ['r','b','c','m','y','k'];
color = [color, color, color, color, color, color, color, color, color, color, color, color];
vert = zeros(0,3);
for j = 1:length(r)
    hold on
%     %uncomment below to also plot ray lines with the shape points
%     for i = 1:10:length(limb_starts{j})
%         line([limb_starts{j}(i,1) limb_ends{j}(i,1)], [limb_starts{j}(i,2) limb_ends{j}(i,2)],[limb_starts{j}(i,3) limb_ends{j}(i,3)],'Color', [1,0,0,0.2],'LineWidth',1);
%     end
    scatter3(edge_points_bc{j}(:,1),edge_points_bc{j}(:,2),edge_points_bc{j}(:,3),'filled')
    vert = [vert;(edge_points_bc{j})];
    drawnow;
    %disp(j)
    hold off    
     
end
axis('equal')
view(45,45)

%%
[~, shapePnts, shapePntNhats] = shape_from_limbs_lo(limb_starts, limb_ends, r, 10,1);
toc
%plots the shapePnt results of shape_from_limbs
ptCloud = pointCloud(shapePnts,'Normal',shapePntNhats);
figure(6)
pcshow(ptCloud)

xlabel('X Axis')

ylabel('Y Axis')
zlabel('Z Axis')

%% remove outliers
%downsample
ptCloud_downsample = pcdownsample(ptCloud,'gridAverage',0.001);

%filter out noisy particles
[ptCloud2,inlierIndices,outlierIndices] = pcdenoise(ptCloud_downsample,'NumNeighbors',10,'Threshold',5);

figure(7)
pcshow(ptCloud2)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
disp('Outliers: ')
disp(length(outlierIndices))

%% 
p = [ptCloud2.Location(:,1), ptCloud2.Location(:,2), ptCloud2.Location(:,3)];
% p = [pnts(:,1),pnts(:,2),pnts(:,3)];
[t,tnorm]=MyRobustCrust(p);

figure()
title('Output Triangulation','fontsize',14)
axis equal
trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b');%'facecolor',[0.0 0.0 0.1,0.1],'edgecolor',[0.5 0.5 0.5 0.1]);%plot della superficie trattata
view(3);
axis('equal')
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')

%% save shape model to PLY file
h = trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b');axis equal
plywrite(datasetPath+body+"_model.ply",h.Faces,h.Vertices)
