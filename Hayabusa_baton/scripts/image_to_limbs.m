%%
% 
%
% 
%takes image data and puts all rays in first image's BC coordinate frame
%Dahlia Baker
%Last edit - December 2 2019
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

function [limb_starts, limb_ends, test_long] = image_to_limbs(img_list, mask_list,  lat_list, long_list, z_list, fov_angle)
  
    j = 1;
    
    while j <= length(img_list)

        asteroid = imread(img_list(j));
        asteroid = rgb2gray(asteroid);
        masksize = mask_list(j,:);

        %hardcoded for every iteration
        sample_num = 25;
        itr_num = 1500;
        angdiff = 5;
        
        [trim_u, trim_v] = edge_finding(asteroid,masksize, sample_num, itr_num, angdiff);

        z = z_list(j);

        [edge_points_bc{j}, ~] = edge_to_3d(z, fov_angle, trim_u, trim_v);

        lat = lat_list(j);%degrees
        long = long_list(j);%degrees
        
        
        %did I confuse lat and long?
        %if j ~= 1
            theta = -deg2rad(long);%-deg2rad(long_list(1));
            phi   = deg2rad(lat);%-deg2rad(lat_list(1));
        
            T_y = [cos(theta) 0 -sin(theta);...
                   0 1 0;...
                   sin(theta) 0 cos(theta)];
           
%             T_z = [cos(phi) sin(phi) 0;...
%                    -sin(phi) cos(phi) 0;...
%                    0 0 1];
            T_x = [1 0 0;...
                   0 cos(phi) -sin(phi);...
                   0 sin(phi) cos(phi)];
            
            new_pts = zeros(3,length(edge_points_bc{j}(:,1)));

            for i = 1:length(edge_points_bc{j}(:,1))
                new_pts(:,i) = (T_x*T_y)*(edge_points_bc{j}(i,:)');
            end

            new_pts = new_pts';

            rotated_z = ((T_x*T_y)*[0;0;.5])';
            clear new_rays_bc

        %elseif j == 1
            
%             rotated_z = [0; 0; .5]';
%             new_pts = edge_points_bc{j};
%             theta = 0;
%         end

        new_rays_bc(:,1:3) = new_pts(:,1:3) + rotated_z;
        new_rays_bc(:,4:6) = new_pts(:,1:3) - rotated_z;


%         len = length(new_rays_bc(:,1));

%         for i = 1:len-1
%             line([new_rays_bc(i,1) new_rays_bc(i,4)], [new_rays_bc(i,2) new_rays_bc(i,5)], [new_rays_bc(i,3) new_rays_bc(i,6)],'Color',color(j))
%         end

    
        limb_starts{j} = new_rays_bc(:,1:3);
        limb_ends{j}  = new_rays_bc(:,4:6);
        
        
        for i = 1:length(limb_starts{j}(:,1))
            limb_starts{j}(i,:) = [limb_starts{j}(i,3),-limb_starts{j}(i,1), -limb_starts{j}(i,2)];
            limb_ends{j}(i,:) = [limb_ends{j}(i,3), -limb_ends{j}(i,1), -limb_ends{j}(i,2)];
        end
    
        test_long(j)   = theta;
        
        disp('image number processed: ')
        disp(j)
        j = j+1;
        
    end

end