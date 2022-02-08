%Run Script - OSIRIS-REx Version

% also very done with this
%Dahlia Baker
%Last Modified: 10/18/2021

clear all
close all
clc
tic
%set initial values
load_script = 'orex/orex_test_data.mat';
img_path = ["orex/To_Send"];
fov_angle = 0.8;
ext = 5; %length of ray extension
limb = 1; %1 for include terminator, 0 for limb-only
IR = 0;%0 if no IR data, 1 if yes
density = 5; %number of degrees between silhouette sampled points
body = 'bennu';
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
fov_angle = 0.8;

for i = 1:length(img_path)
    addpath(img_path(i));   
end

% for i = 1:69
%     CB(:,:,i) = rotz(180)*CB(:,:,i);
% end

%%

[limb_starts, limb_ends, edge_points_bc] = image_to_limbs_orex(names, r, fov_angle,CB,sun_pos,phase,limb,ext,IR,IRdat,density,body,ir_list);

%%
figure(5)
view([45 30])
hold on
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')
color = ['r','b','c','m','y','k'];
color = [color, color, color, color, color, color, color, color, color, color, color, color];

for j = 1:69
    %uncomment below to also plot ray lines with the shape points
     for i = 1:72
        
         line([limb_starts{j}(i,1) limb_ends{j}(i,1)], [limb_starts{j}(i,2) limb_ends{j}(i,2)],[limb_starts{j}(i,3) limb_ends{j}(i,3)],'Color', [1,0,0,0.2],'LineWidth',1);
     end
%      drawnow;
%         for i = 1:108
%             line([limb_starts{j}(i,1) limb_ends{j}(i,1)], [limb_starts{j}(i,2) limb_ends{j}(i,2)],[limb_starts{j}(i,3) limb_ends{j}(i,3)],'Color', color(j),'LineWidth',1);
%         end
  scatter3(edge_points_bc{j}(1:36,1),edge_points_bc{j}(1:36,2),edge_points_bc{j}(1:36,3),'filled')
  drawnow;  
  disp(j)
  scatter3(edge_points_bc{j}(37:72,1),edge_points_bc{j}(37:72,2),edge_points_bc{j}(37:72,3),'filled')
  drawnow;  
  disp(j)
%   
        
end
axis('equal')

%%
clear limb_starts_sub limb_ends_sub l_s l_e
sub = 1:69;
k = 1;
%figure(5);

for i = 1:69
    limb_starts_sub{k} = limb_starts{i};
    limb_ends_sub{k} = limb_ends{i};
    
%     for j = 1:72
%         
%          line([limb_starts_sub{k}(j,1) limb_ends_sub{k}(j,1)], [limb_starts_sub{k}(j,2) limb_ends_sub{k}(j,2)],[limb_starts_sub{k}(j,3) limb_ends_sub{k}(j,3)],'Color', [1,0,0,0.2],'LineWidth',1);
%     end
%      drawnow;
     k = k+1;
end
hold off
l_s = flip(limb_starts_sub);
l_e = flip(limb_ends_sub);
%%

%figure(7)
hold on
%k = 1;
for k = [1:10]
     for j = 1:72
        if k == 1 
         line([l_s{k}(j,1) l_e{k}(j,1)], [l_s{k}(j,2) l_e{k}(j,2)],[l_s{k}(j,3) l_e{k}(j,3)],'Color', [0,0,1,0.2],'LineWidth',1);
        else
            line([l_s{k}(j,1) l_e{k}(j,1)], [l_s{k}(j,2) l_e{k}(j,2)],[l_s{k}(j,3) l_e{k}(j,3)],'Color', [1,0,0,0.2],'LineWidth',1);
        end
         
     end
     scatter3(l_s{k}(:,1),l_s{k}(:,2),l_s{k}(:,3),'g')
     scatter3(l_e{k}(:,1),l_e{k}(:,2),l_e{k}(:,3),'b')
    drawnow;
    disp(k)
    %k = k+1;
end
%%
[shapeEndPnts,shapePnts, shapePntNhats] = shape_from_limbs_orex(limb_starts,limb_ends, r, 10,2,0,[0,0,-1]);
toc
%plots the shapePnt results of shape_from_limbs
ptCloud = pointCloud(shapePnts,'Normal',shapePntNhats);
figure(8)
pcshow(ptCloud)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')

%%
%filter out noisy particles
[ptCloud,inlierIndices,outlierIndices] = pcdenoise(ptCloud2,'NumNeighbors',10,'Threshold',5);

figure(5)
pcshow(ptCloud3)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
disp('Outliers: ')
disp(length(outlierIndices))

%%
p = [ptCloud.Location(:,1), ptCloud.Location(:,2), ptCloud.Location(:,3)];

[t,tnorm]=MyRobustCrust(p);

figure(6)
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
plywrite('orex/orex.ply',h.Faces,h.Vertices);

%%
%add 'p'  and 't' and 'tnorm'  back into line below when you have a shape model to save
%save('0phase_sim_bennu/0deg_bennu_output_data.mat','limb_ends','limb_starts','ptCloud','shapeEndPnts','shapePntNhats','shapePnts','t','tnorm','p','-v7.3')



