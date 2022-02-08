%%
% Update
% 2/3/2022
% Dahlia Baker
% limb only processing




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

function [limb_starts, limb_ends, edge_points_bc] = image_to_limbs_lo(img_list, z_list, fov_angle,CB,sun_pos,phase,ext)
  
    j = 1;
    
    while j <= length(img_list)

        asteroid = imread(img_list(j));
        asteroid = rgb2gray(asteroid); 
        asteroid(asteroid<uint8(10)) = uint8(0);
        asteroid = asteroid*1000;
          
        [trim_u, trim_v,~,~,mid_pt_u,mid_pt_v] = edge_finding_lo(asteroid);
        
        %check sign of y comp of SunB
        cam_pos = [0,0,z_list(j)]; %uncomment for regular cases
        sunb = cam_pos - (CB(:,:,j)'*sun_pos(j,:)')';
        %sunb = cam_pos(j,:) - (CB(:,:,j)'*sun_pos(j,:)')';
        if sunb(2) >= 0
            dir = 1;
        else
            dir = -1;
        end
        
        [edge_points{j}, edge_rays{j}, new_trim_u,new_trim_v] = edge_to_3d_lo(z_list(j), fov_angle, trim_u, trim_v,sunb,mid_pt_u,mid_pt_v,dir,ext);
        %plot them one over another
       
        if j ==0
            %ast_flip = flip(asteroid,1);
            figure(1)
            imshow(asteroid)
            hold on
            grid on
            scatter(new_trim_u,new_trim_v,'filled','b')
            scatter(mid_pt_v,mid_pt_v,'filled','r')
            %plot lines from sun direction
            %asteroid = imcrop(asteroid);
            legend({'limb','center'},'FontSize',24)
            xlabel('X (pixels)','FontSize',16)
            ylabel('Y (pixels)','FontSize',16)
            title(string(j),'FontSize',24)
            hold off    
        end
        
        T = CB(:,:,j);
   
        [edge_points_bc{j}, new_rays_bc] = CFtoBF_orex(T, edge_points{j}, edge_rays{j});

        limb_starts{j} = new_rays_bc(:,1:3);
        limb_ends{j}  = new_rays_bc(:,4:6);
       
        j = j+1;
    end
end