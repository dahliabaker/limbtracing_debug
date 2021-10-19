%%
% 
%New Image_to_limbs code without box-shrinking
% 
%takes image data and puts all rays in first image's BC coordinate frame
%Dahlia Baker
%Last edit - January 6 2020
%
% Input:
% img_list - 1 x n list of string names of images to be processed
%
% mask_list - 
%       n x 4 array of mask parameters in format [vmin vmax umin
%       umax]. Has to be defined for each image or it wont converge
%
% lat_list - 1 x n list of latitudes in center of images
%
% long_list - 1 x n list of longitudes in center of images
%
% z_list - 1 x n list of z distances (camera to body) in km
%
% fov_angle - a camera constant, how many km each pixel represents at
% bodycentric origin
%
%

function [limb_starts, limb_ends, edge_points_bc] = image_to_limbs_orex(img_list, z_list, fov_angle,CB,sun_pos,phase,limb,ext,IR,IRdat,density,body,ir_list)
  
    j = 1;
    
    while j <= length(img_list)

        asteroid = imread(img_list(j));
        
        %options for simulated images
         asteroid = rgb2gray(asteroid); %toggle on or off based on input
         if j == 49
             asteroid(asteroid<uint8(20)) = uint8(0);
         else
            asteroid(asteroid<uint8(10)) = uint8(0);
         end
         asteroid = asteroid*1000;
     
        [trim_u, trim_v,~,~,mid_pt_u,mid_pt_v] = edge_finding_canny(asteroid, density,limb);
        
        %check sign of y comp of SunB
        cam_pos = [0,0,z_list(j)]; %uncomment for regular cases
        sunb = cam_pos - (CB(:,:,j)'*sun_pos(j,:)')';
        %sunb = cam_pos(j,:) - (CB(:,:,j)'*sun_pos(j,:)')';
        if sunb(2) >= 0
            dir = 1;
        else
            dir = -1;
        end
        vis_limb = 0;
        [edge_points{j}, edge_points_t{j}, edge_rays{j}, edge_rays_t{j}, new_trim_u,new_trim_v, new_term_u,new_term_v] = edge_to_3d_orex(z_list(j), fov_angle, trim_u, trim_v,sunb,mid_pt_u,mid_pt_v,dir,phase(j),limb,ext,vis_limb);
        %plot them one over another
       

        if j >0
            %ast_flip = flip(asteroid,1);
            figure(1)
            imshow(asteroid)
            hold on
            grid on
            scatter(new_trim_u,new_trim_v,'filled','b')
            scatter(new_term_u,new_term_v,'filled','r')
            %plot lines from sun direction
            %asteroid = imcrop(asteroid);
            legend({'limb','terminator'},'FontSize',24)
            xlabel('X (pixels)','FontSize',16)
            ylabel('Y (pixels)','FontSize',16)
            title(string(j),'FontSize',24)
            hold off
            
        end
        
        T = CB(:,:,j);
    %commented out for Hayabusa case
        sun_v3 = sun_pos(j,:)+[0,0,-100];%in camera frame
        rhat = [0,0,-1];%camera pointing vector in camera frame
        theta_yz = asind(norm(cross(sun_v3,rhat))/(norm(sun_v3)*norm(rhat))); %angle between rhat and sun vector in y-z plane
        %check for which side of camera sun is on
        if sun_v3(2) >= 0
            ang = -theta_yz;
        else 
            ang = theta_yz;
        end
        rot1 = [1 0 0; 0 cosd(ang) sind(ang); 0 -sind(ang) cosd(ang)];%ang rotation about first axis
        theta_xz = 180-atan2d(sun_v3(1),sun_v3(3));
        rot2 = [cosd(theta_xz) 0 sind(theta_xz); 0 1 0; -sind(theta_xz) 0 cosd(theta_xz)]; 
        
        rot = rot2*rot1;
        
        [edge_points_bc{j}, new_rays_bc] = CFtoBF_orex(T, rot, edge_points{j}, edge_points_t{j}, edge_rays{j},edge_rays_t{j});


        limb_starts{j} = new_rays_bc(:,1:3);
        limb_ends{j}  = new_rays_bc(:,4:6);
%         limb_starts{j} = new_rays_bc(:,4:6);
%         limb_ends{j}  = new_rays_bc(:,1:3);


        j = j+1;
    end
