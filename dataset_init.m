%dataset_init.m
%10/7/2021

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
    r(x) = C(1); %distance from body to camera
end

%change save file name to match case and path
save('72phase_sim_bennu/bennu_72.mat','phase','CB','sun_pos','cam_pos','img_name','r');
