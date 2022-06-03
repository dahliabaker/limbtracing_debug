%Run Script
%2/16/2022
%Dahlia Baker
clear all
close all
clc

restoredefaultpath
%%
%set initial values
addpath('scripts_limbonly')
datasetPath = "new_datasets/45bennu_sim/";
img_path = datasetPath+"bennu_automated_images";%;"mv_bennu/bennu_automated_images_2"];
fov_angle = 0.8;
ext = 1; %length of ray extension
limb = 0; %1 for include terminator, 0 for limb-only
load(datasetPath+"bennu_72.mat");
    
for i = 1:length(img_name)
    img_name(i) = img_path+"/render"+string(i)+".png";
end

clear cam_pos
[limb_starts, limb_ends, edge_points_bc] = image_to_limbs_term(img_name, r, fov_angle,CB,sun_pos,ext);

%%
figure(4)
hold on
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')
color = ['r','b','c','m','y','k'];
color = [color, color, color, color, color, color, color, color, color, color, color, color];

for j = 1:72
    for i = 26:50
        line([limb_starts{j}(i,1) limb_ends{j}(i,1)], [limb_starts{j}(i,2) limb_ends{j}(i,2)],[limb_starts{j}(i,3) limb_ends{j}(i,3)],'Color', [1,0,0,0.2],'LineWidth',1);
    end
        scatter3(edge_points_bc{j}(26:50,1),edge_points_bc{j}(26:50,2),edge_points_bc{j}(26:50 ,3),'filled')
        drawnow; 
%         disp('stop')
        
end
axis('equal')
%%
[~,shapePnts, shapePntNhats] = shape_from_limbs_terminator(limb_ends,limb_starts,r, 10,25);

%plots the shapePnt results of shape_from_limbs
ptCloud = pointCloud(shapePnts,'Normal',shapePntNhats);
figure(6)
pcshow(ptCloud,'MarkerSize',1)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')

%%
hold on
color = ['r','b','c','m','y','k'];
color = [color, color, color, color, color, color, color, color, color, color, color, color];

for j = 1:10
%     for i = 1:20
%         line([limb_starts{j}(i,1) limb_ends{j}(i,1)], [limb_starts{j}(i,2) limb_ends{j}(i,2)],[limb_starts{j}(i,3) limb_ends{j}(i,3)],'Color', [1,0,0,0.2],'LineWidth',1);
%     end
        scatter3(edge_points_bc{j}(:,1),edge_points_bc{j}(:,2),edge_points_bc{j}(:,3),8,'r','filled')
        drawnow; 
%         disp('stop')
%         
end
axis('equal')

%%
%downsample
%ptCloud2 = pcdownsample(ptCloud2,'gridAverage',0.005);

%filter out noisy particles
[ptCloud3,inlierIndices,outlierIndices] = pcdenoise(ptCloud,'NumNeighbors',10,'Threshold',2);

figure(7)
pcshow(ptCloud3)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
disp('Outliers: ')
disp(length(outlierIndices))

%%

p = [ptCloud3.Location(:,1), ptCloud3.Location(:,2), ptCloud3.Location(:,3)];
t = boundary(p(:,1),p(:,2),p(:,3));  
h = trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b');

%%
%plywrite("45itokawa_model.ply",h.Faces,h.Vertices);
