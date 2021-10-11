%Run Script
%6/29/2021
%Dahlia Baker
clear all
close all
clc
tic
%set initial values
load_script = 'orex/orex_test_data.mat';
img_path = ["orex/To_Send"];%;"mv_bennu/bennu_automated_images_2"];
fov_angle = 0.8;
%phase = 10; %phase of test case - can be an
ext = 1; %length of ray extension
limb = 1; %1 for include terminator, 0 for limb-only
IR = 0;%0 if no IR data, 1 if yes
density = 5; %number of degrees between silhouette sampled points
body = 'itoka';
if IR == 0
    IRdat = [];
else
    load('10phase_sim_itokawa/IRedgedata.mat')
    IRdat = irSilhoutte;
end
    
%%
load(load_script);
addpath('scripts/');

%load('rotation_sub_10.mat');
%load('rotation_sub.mat')
%load('rotation_3deg.mat')

for i = 1:length(img_path)
    addpath(img_path(i));   
end
    
[limb_starts, limb_ends, edge_points_bc] = image_to_limbs_orex(names, r, fov_angle,CB,sun_pos,phase,limb,ext,IR,IRdat,density,body);

figure(2)
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
%%
[shapeEndPnts,shapePnts, shapePntNhats] = shape_from_limbs_orex(limb_starts,limb_ends, r, 10,2,0,[0,0,-1]);
toc
%plots the shapePnt results of shape_from_limbs
ptCloud = pointCloud(shapePnts,'Normal',shapePntNhats);
figure(3)
pcshow(ptCloud)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')

%%
shapePntsNew = [];
shapePntNhatsNew = [];
spin_pole = [0,0,-1];

for i = 1:length(shapePnts(:,1))
    nhat = shapePntNhats(i,:);
    vec_diff = nhat-spin_pole;

    if vec_diff(3) < 1e-2 || vec_diff(3) > (2-1e-2)
        %skip this ray, too close to spin pole
        %remove these two lines later, just taking out the pole
        %adjustment functionality to test
    else
        shapePntsNew = [shapePntsNew; shapePnts(i,:)];
        shapePntNhatsNew = [shapePntNhatsNew; shapePntNhats(i,:)];
    end
end

figure
h = histogram(shapePntsNew(:,3),100);
title('Point Z Spread','FontSize', 18)
xlabel('Z coord (km)','FontSize', 16)
ylabel('Count','FontSize', 16)

shapePntsNewNew = [];
shapePntNhatsNewNew = [];
shapePntsNewNewNew =[];
shapePntNhatsNewNewNew = [];

out = find(h.Values == 0);
half = median(h.Values);
if length(out)>1
    data_bin = round(length(shapePntsNew(:,3))/100);
   
    first = min(out)*data_bin*0.5;
    second = out(end)*data_bin;
    % for i = 1:length(out)
    %     if out(i) > half
    %         disp(out(i))
    %         second = out(i)*data_bin;
    %         return
    %     end
    % end

%find a skip in the histogram, and don't start saving until after the
%first and before the last
%sort points by Z

    [~,I] = sort(shapePntsNew(:,3));

    for i = 1:length(shapePntsNew(:,3))
        shapePntsNewNew = [shapePntsNewNew; shapePntsNew(I(i),:)];
        shapePntNhatsNewNew = [shapePntNhatsNewNew; shapePntNhatsNew(I(i),:)];
    end
%find all data below out(1) && less than median


    shapePntsNewNewNew = shapePntsNewNew(first:second,:);
    shapePntNhatsNewNewNew = shapePntNhatsNewNew(first:second,:);
else
    shapePntsNewNewNew = shapePntsNew;
    shapePntNhatsNewNewNew = shapePntNhatsNew;
end


%%
ptCloud = pointCloud(shapePntsNewNewNew,'Normal',shapePntNhatsNewNewNew);
figure(4)
pcshow(ptCloud2)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')



%%
%downsample
%ptCloud2 = pcdownsample(ptCloud2,'gridAverage',0.005);

%filter out noisy particles
[ptCloud3,inlierIndices,outlierIndices] = pcdenoise(ptCloud,'NumNeighbors',10,'Threshold',2);

figure(5)
pcshow(ptCloud3)
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
disp('Outliers: ')
disp(length(outlierIndices))

%%
% k = 1;
% for i = 1:72
%     for j = 1:36
%         pnts(k,1:3) = [edge_points_bc{i}(j,1);edge_points_bc{i}(j,2);edge_points_bc{i}(j,3)];
%         k = k+1;
%     end
% end
% points = [];
p = [ptCloud.Location(:,1), ptCloud.Location(:,2), ptCloud.Location(:,3)];
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

%%
h = trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b');
plywrite('10phase_sim_itokawa/10itokawa_IR.ply',h.Faces,h.Vertices);

%%
%add 'p'  and 't' and 'tnorm'  back into line below when you have a shape model to save
%save('0phase_sim_bennu/0deg_bennu_output_data.mat','limb_ends','limb_starts','ptCloud','shapeEndPnts','shapePntNhats','shapePnts','t','tnorm','p','-v7.3')



