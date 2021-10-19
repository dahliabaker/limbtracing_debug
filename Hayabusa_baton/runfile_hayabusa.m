%Run Script - Hayabusa Version
% i am done with this

%10/18/2021
%Dahlia Baker
clear all
close all
clc
tic
%set initial values
load_script = 'hayabusa/hayabusa.mat';
img_path = ["hayabusa/hayabusa_images"];
fov_angle = 0.8;
ext = 1; %length of ray extension
limb = 1; %1 for include terminator, 0 for limb-only
IR = 0;%0 if no IR data, 1 if yes
density = 5; %number of degrees between silhouette sampled points
body = 'itoka';
if IR == 0
    IRdat = [];
    ir_list = [];
else
    load('10phase_sim_itokawa/IRedgedata.mat')
    IRdat = irSilhoutte;
    load('ir_list_i.mat')
    ir_list = ir_imgs;
end
    
%%
load(load_script);
addpath('scripts/');

for i = 1:length(img_path)
    addpath(img_path(i));   
end


%hayabusa image set conditioning
good = 20:1:75;


[limb_starts, limb_ends, edge_points_bc] = image_to_limbs_orex(img_name(good), r(good), fov_angle,CB(:,:,good),sun_pos(good,:),phase(good),limb,ext,IR,IRdat,density,body,ir_list);

%%
figure(3)
hold on
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')
color = ['r','b','c','m','y','k'];
color = [color, color, color, color, color, color, color, color, color, color, color, color];

for j = 1:length(good)
    %uncomment below to also plot ray lines with the shape points
%      for i = 1:4:72
%          line([limb_starts{j}(i,1) limb_ends{j}(i,1)], [limb_starts{j}(i,2) limb_ends{j}(i,2)],[limb_starts{j}(i,3) limb_ends{j}(i,3)],'Color', [1,0,0,0.2],'LineWidth',1);
%      end
%         for i = 1:108
%             line([limb_starts{j}(i,1) limb_ends{j}(i,1)], [limb_starts{j}(i,2) limb_ends{j}(i,2)],[limb_starts{j}(i,3) limb_ends{j}(i,3)],'Color', color(j),'LineWidth',1);
%         end
        scatter3(edge_points_bc{j}(:,1),edge_points_bc{j}(:,2),edge_points_bc{j}(:,3),'filled')
        drawnow;
        %disp(j)  
end
axis('equal')


%%
[shapeEndPnts,shapePnts, shapePntNhats] = shape_from_limbs_orex(limb_starts,limb_ends, r(good), 10,5,0,[0,0,-1]);
toc
%plots the shapePnt results of shape_from_limbs
ptCloud = pointCloud(shapePnts,'Normal',shapePntNhats);
figure(4)
pcshow(ptCloud)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')



%%
%filter out noisy particles
[ptCloud,inlierIndices,outlierIndices] = pcdenoise(ptCloud,'NumNeighbors',10,'Threshold',2);

figure(5)
pcshow(ptCloud)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
disp('Outliers: ')
disp(length(outlierIndices))

%%
p = [ptCloud.Location(:,1), ptCloud.Location(:,2), ptCloud.Location(:,3)];

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

%%
h = trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b');
plywrite('10phase_sim_itokawa/10itokawa_IR.ply',h.Faces,h.Vertices);

%%
%add 'p'  and 't' and 'tnorm'  back into line below when you have a shape model to save
%save('0phase_sim_bennu/0deg_bennu_output_data.mat','limb_ends','limb_starts','ptCloud','shapeEndPnts','shapePntNhats','shapePnts','t','tnorm','p','-v7.3')



