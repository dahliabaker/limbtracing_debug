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

function [limb_starts, limb_ends, edge_points_bc] = image_to_limbs_canny(img_list, lat_list, long_list, z_list, fov_angle)
  
    j = 1;
    
    while j <= length(img_list)

        asteroid = imread(img_list(j));
        asteroid = rgb2gray(asteroid); %toggle on or off based on input
        asteroid = imadjust(asteroid,[0 1]);
        asteroid = imgaussfilt(asteroid,4);
        %imwrite(asteroid,'blur_'+string(j)+'.png')
        
        angdiff = 2; 
        
        [trim_u, trim_v,E_u,E_v] = edge_finding_canny(asteroid, angdiff);
        
%         figure()
%         gcf = scatter(trim_u,trim_v,'filled');
%         xlabel('x (km)')
%         ylabel('y (km)')
%         xlim([-500,500])
%         ylim([-300,300])
%         saveas(gcf,'edge_'+string(j)+'.png');
        
        z = z_list(j);

        [edge_points{j}, edge_rays{j}] = edge_to_3d(z, fov_angle, trim_u, trim_v);

        lat = lat_list(j);%degrees
        long = long_list(j);%degrees
        
        [edge_points_bc{j}, new_rays_bc] = CFtoBF(lat, long, edge_points{j}, edge_rays{j});
        
        



        limb_starts{j} = new_rays_bc(:,1:3);
        limb_ends{j}  = new_rays_bc(:,4:6);
        
        
        disp('image number processed: ')
        disp(j)
        
        
        j = j+1;
     end
    

end