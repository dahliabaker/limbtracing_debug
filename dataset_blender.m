%Creating data sets
% To go with automated blender images

%camera frame is stationary
%but to body, camera is always rotating, as well as sun

%number of images
n = 72;

%camera vector in initial body frame
C = [100,0,0];

%sun vector in initial body frame (equivalent to values in blender)
S = [200,0,0];

%amount of angular change (degree) btwn images
rotate_by = 5;
start_angle = -5;

%calculate and spit out the phase angle
phase = acosd(dot(C,S)/(norm(C)*norm(S)));
disp('phase angle = ')
disp(phase)

%%%initializing the camera frame weirdness
%-90 rotation about second axis
rot1 = [cosd(-90) 0 sind(-90); 0 1 0; -sind(-90) 0 cosd(-90)];
% +/- 180 rotation about third axis
rot2 = [cosd(180), sind(180), 0; -sind(180), cosd(180), 0; 0, 0, 1];
%account for y being in the top corner
rot3 = [1 0 0; 0 -1 0; 0 0 1];

for x = 1:n
    %current angle of rotation
    angle = (start_angle * (pi/180)) + (x*1) * (rotate_by * (pi/180)); %degree to radian
    eul = [angle,0,0]; %euler rotation
    rotmZYX(:,:,x) = eul2rotm(eul); %the rotation matrix applied to the body
    
    %body-to-camera frame rotation
    CB(:,:,x) = (rot3*rot2*rot1*(rotmZYX(:,:,x)));
    sun_pos(x,:) = [-S(3)+C(3),-S(2)+C(2),-S(1)+C(1)];%sun position in camera frame
    cam_pos(x,:) = [0,0,0];%camera_position in camera frame
    %change img_name path to point towards wherever you saved your images
    img_name(x) = '72phase_sim_bennu/bennu_automated_images/render'+string(x)+'.png';
    
    r(x) = C(1);
end


%%
%camera frame is stationary
%but to body, camera is always rotating, as well as sun
n = 360; %number of images
%initially
C = [100,0,0]; %camera vector in initial body frame
%S = [-200,0,0];
%S = [200,0,0];
S = [190,33.502,0];
%S = [165,119.8793,0]; %sun vector in initial body frame, phase = 36.8699
%S = [200,35.2655, 0];%phase = 8.9270
%S = [195,-60,60];
rotate_by = 5;
start_angle = -5;    
rotmZYX = zeros(3,3,n);

%calculate and spit out the phase angle
phase = acosd(dot(C,S)/(norm(C)*norm(S)));


%%%initializing the camera frame weirdness
%-90 rotation about second axis
rot1 = [cosd(-90) 0 sind(-90); 0 1 0; -sind(-90) 0 cosd(-90)];
% +/- 180 rotation about third axis
rot2 = [cosd(180), sind(180), 0; -sind(180), cosd(180), 0; 0, 0, 1];
%account for y being in the top corner
rot3 = [1 0 0; 0 -1 0; 0 0 1];

for x = 1:n
    %current angle of rotation
    angle = (start_angle * (pi/180)) + (x*1) * (rotate_by * (pi/180)); %degree to radian
    eul = [angle,0,0]; %euler rotation
    rotmZYX(:,:,x) = eul2rotm(eul); %the rotation matrix applied to the body
    
    %body-to-camera frame rotation
    CB(:,:,x) = (rot3*rot2*rot1*(rotmZYX(:,:,x)));
    sun_pos(x,:) = [-S(3)+C(3),-S(2)+C(2),-S(1)+C(1)];%sun position in camera frame
    cam_pos(x,:) = [0,0,0];%camera_position in camera frame
    img_name(x) = '72phase_sim_bennu/bennu_automated_images/render'+string(x)+'.png';
    
    r(x) = C(1);
end

%% SOMETHING IS WRONG WITH THIS
%negative euler angle rotation about third axis is the first rotation
% %+90 rotation about first axis
% rot1 = [1 0 0; 0 cosd(90) sind(90); 0 -sind(90) cosd(90)];
% %+90 rotation about second axis
% rot2 = [cosd(180) 0 sind(180); 0 1 0; -sind(180) 0 cosd(180)];
% %-90 rotation about third axis
% rot3 = [cosd(270), sind(270), 0; -sind(270), cosd(270), 0; 0, 0, 1];


%-90 rotation about second axis
rot1 = [cosd(-90) 0 sind(-90); 0 1 0; -sind(-90) 0 cosd(-90)];
% +/- 180 rotation about third axis
rot2 = [cosd(180), sind(180), 0; -sind(180), cosd(180), 0; 0, 0, 1];
%account for y being in the top corner
rot3 = [1 0 0; 0 -1 0; 0 0 1];
%rot3 = eye(3,3);
%%
n = 72;
C = [100,0,0];
for x = 1:72
    %S = phase_angles(:,x);
    angle = (start_angle * (pi/180)) + (x*1) * (rotate_by * (pi/180)); %degree to radian
    eul = [angle,0,0]; %euler rotation
    rotmZYX(:,:,x) = eul2rotm(eul); %the rotation matrix applied to the body
    %error(:,:,x) = eul2rotm([1,0,0]);
    %body-to-camera frame rotation
    %CB(:,:,x) = (rot3*rot2*rot1*(error(:,:,x)*rotmZYX(:,:,x))); 
    CB(:,:,x) = (rot3*rot2*rot1*(rotmZYX(:,:,x)));
%     CB_b(:,:,x) = flip(CB_b(:,:,x),1);
%     CB_b(:,:,x) = flip(CB_b(:,:,x),2);
    %sun_pos(x,:) = [-S(3)+C(3),-S(2)+C(2),-S(1)+C(1)];%sun position in camera frame
    %cam_pos(x,:) = [0,0,0];%camera_position in camera frame
    %img_name(x) = '72phase_sim_bennu/bennu_automated_images/render'+string(x)+'.png';
    %img_name(x) = 'psyche_automated_images/render'+string(x)+'.png';
    %r(x) = 100;
end

%save('bennu_data_36.mat','CB','sun_pos','img_name','r')
%save('72phase_sim_bennu/bennu_72.mat','CB','sun_pos','img_name','r','C');
%save('3deg_itokawa.mat','CB','sun_pos','img_name','r')
save('mv_bennu.mat','CB','sun_pos','img_name','r');






%% testing phase angle center assumptions
clear all
close all
clc

addpath('phase_angle_test')
addpath('scripts')
load('sphere_angles.mat')
fov_angle = 0.8;
z = 100;
%load image and run canny
j = 1;

while j <= length(img_name)

    asteroid = imread(img_name(j));
    asteroid = rgb2gray(asteroid); %toggle on or off based on input
    asteroid = imadjust(asteroid,[0 1]);
    %asteroid = imgaussfilt(asteroid,4);
    %imwrite(asteroid,'blur_'+string(j)+'.png')

    S = phase_angles(:,j);
    C = [100,0,0];
    phase = acosd(dot(C,S)/(norm(C)*norm(S)));

   
    [trim_u, trim_v,E_u,E_v,mid_pt_u,mid_pt_v] = edge_finding_canny(asteroid, 4);

    [edge_points{j}, edge_points_t{j}, edge_rays{j}, edge_rays_t{j}, new_trim_u,new_trim_v, new_term_u,new_term_v] = edge_to_3d_orex(z, fov_angle, trim_u, trim_v,sun_pos(j,:),mid_pt_u,mid_pt_v,1,phase);
    %plot them one over another
    offset_u(j) = 512-mid_pt_u;
    offset_v(j) = 512-mid_pt_v;

    %ast_flip = flip(asteroid,1);
    figure()
    imshow(asteroid)
    hold on
    scatter(new_trim_u,new_trim_v,'filled','b')
    scatter(new_term_u,new_term_v,'filled','r')
    %plot lines from sun direction

    legend({'limb','terminator'},'FontSize',12)
    hold off
%     saveas(gcf,'phase_angle_test/output_sphere/edge_'+string(j)+'.png','png')
    if phase >= 90
        disp('stop')
    end
        
    j = j+1;
end

%%
figure()
subplot(2,1,2)
scatter(0:5:160,offset_v(1:33),'filled','b')
ylabel('Y Center Offset (pixels)')
xlabel('Phase Angle (degrees)')
grid on
subplot(2,1,1)
scatter(0:5:160,offset_u(1:33),'filled','r')
ylabel('X Center Offset (pixels)')
grid on
title('Center Offset from Truth versus Phase Angle: Sphere')

%%
phase = 0:5:160;
corrcoef(phase,offset_v(1:33));
[p,errorEst] = polyfit(phase,offset_v(1:33),4);
pop_fit = polyval(p,phase,errorEst);

figure()
plot(phase,pop_fit,'-',phase,offset_v(1:33),'+')
title('Center Offset versus Phase Angle: Sphere')
xlabel('phase angle (degrees)')
ylabel('offset (pixels)')
legend('Polynomial Model','Data','Location','NorthWest');

res = offset_v(1:33)-pop_fit;
figure, plot(phase,res,'+');
title('Residuals for the Quadratic Polynomial Fit')
