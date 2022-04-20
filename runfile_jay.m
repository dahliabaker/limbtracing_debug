%Run Script
%6/29/2021
%Dahlia Baker
clear all
close all
clc
tic
%set initial values
% datasetPath = "new_datasets/15bennu_sim/";img_path = datasetPath+"/bennu_automated_images";body = "bennu";
datasetPath = "new_datasets/15itokawa_sim/";img_path = datasetPath+"/itokawa_automated_images";body = "itokawa";
fov_angle = 0.8;
ext = 1; %length of ray extension (in km)
limbOnly = 0; % 0 to include terminator, 1 for limb-only
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
if limbOnly == 1
    [limb_starts, limb_ends, edge_points_bc_l] = image_to_limbs_lo(img_name,ir_imgs, r, fov_angle,CB,sun_pos,ext);
    termRays_starts = nan;% for plotting
    termRays_ends = nan;% for plotting
    edge_points_bc_t = nan;% for plotting
else
    [limb_starts, ...
     limb_ends, ...
     edge_points_bc_l, ...
     termRays_starts, ...
     termRays_ends, ...
     edge_points_bc_t] = image_to_limbs_term(img_name, r, fov_angle,CB,sun_pos,ext);
end
% 
%% plot sampled pixels and associated rays
figure(5),clf
hold on
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')
color = ['r','b','c','m','y','k'];
color = [color, color, color, color, color, color, color, color, color, color, color, color];
% subplot(122),p = patch('Faces',obj.f.v,'Vertices',obj.v);
p.FaceColor = 'none';
view(45,20)
axis equal
grid on, grid minor, hold on 

for j = 1:length(r)
    % limb rays and point plot
    hold on, grid on, grid minor, title('Limb Points and Rays')
%     imshow(".\new_datasets\15itokawa_sim\itokawa_automated_images\render"+num2str(j)+".png")
    xlabel('X Axis'),ylabel('Y Axis'),zlabel('Z Axis')
    % sampled limb points
    scatter3(edge_points_bc_l{j}(:,1),edge_points_bc_l{j}(:,2),edge_points_bc_l{j}(:,3),'filled','b')
    axis equal, view(45,45)
%     % plot limb rays
%     for i = 1:length(limb_starts{j})
%         line([limb_starts{j}(i,1) limb_ends{j}(i,1)], [limb_starts{j}(i,2) limb_ends{j}(i,2)],[limb_starts{j}(i,3) limb_ends{j}(i,3)],'Color','g','LineWidth',1);
%     end
    % term rays and point plot
%     subplot(122), hold on, grid on, grid minor, title('Terminator Points and Rays')
    xlabel('X Axis'),ylabel('Y Axis'),zlabel('Z Axis')
    if iscell(termRays_ends)
        % sampled terminator points
        scatter3(edge_points_bc_t{j}(:,1),edge_points_bc_t{j}(:,2),edge_points_bc_t{j}(:,3),'filled','r')
%         axis equal, view(90,0)
%         text(mean(edge_points_bc_t{j}(:,1)),mean(edge_points_bc_t{j}(:,2)),mean(edge_points_bc_t{j}(:,3)),num2str(j))
%         % plot term rays
%         for i = 1:length(termRays_starts{j})
% %             line([termRays_starts{j}(i,1) termRays_ends{j}(i,1)], [termRays_starts{j}(i,2) termRays_ends{j}(i,2)],[termRays_starts{j}(i,3) termRays_ends{j}(i,3)],'Color','r','LineWidth',1);
%             line([termRays_starts{j}(i,1) -1/10*termRays_starts{j}(i,1)], [termRays_starts{j}(i,2) -1/10*termRays_starts{j}(i,2)],[termRays_starts{j}(i,3) -1/10*termRays_starts{j}(i,3)],'Color','r','LineWidth',1);
%         end
    end
    
    drawnow
end

%% get shape model points from rays
% set up parallel pool
gcp();

% get shape points
if limbOnly == 1
    shapePnts = shape_from_limbPatches(limb_starts,limb_ends);% limb only processing
else
    shapePnts = shape_from_limb_and_term(limb_starts,...
                                         limb_ends,...
                                         termRays_starts,...
                                         termRays_ends);
end
toc

%plots raw shape point results
ptCloud = pointCloud(shapePnts);
figure(6)
pcshow(ptCloud)
xlabel('X Axis'),ylabel('Y Axis'),zlabel('Z Axis')
%% denoise the model

%filter out noisy particles
[ptCloud2,inlierIndices,outlierIndices] = pcdenoise(ptCloud_downsample,'NumNeighbors',10,'Threshold',5);

figure(7)
pcshow(ptCloud2)
xlabel('X Axis'),ylabel('Y Axis'),zlabel('Z Axis')
disp('Outliers: ')
disp(length(outlierIndices))

%% generate 
p = [ptCloud2.Location(:,1), ptCloud2.Location(:,2), ptCloud2.Location(:,3)];
t = boundary(p(:,1),p(:,2),p(:,3));  

figure()
title('Output Triangulation','fontsize',14)
axis equal
h = trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b');
view(3);
axis('equal')
xlabel('X Axis'),ylabel('Y Axis'),zlabel('Z Axis')

%% save shape model to PLY file
plywrite(datasetPath+body+"_model.ply",h.Faces,h.Vertices)
