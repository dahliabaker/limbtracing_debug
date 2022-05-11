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

function [limb_starts, limb_ends, edge_points_bc] = image_to_limbs_lo(img_list,ir_imgs, z_list, fov_angle,CB,sun_pos,ext)
  
    j = 1;
    
    while j <= length(img_list)
        
        % Optical Image
        asteroid = imread(img_list(j));
        asteroid = rgb2gray(asteroid); 
        asteroid(asteroid<uint8(20)) = uint8(0);
        asteroid = asteroid*1000;
          
        [trim_u, trim_v,~,~,mid_pt_u,mid_pt_v] = edge_finding_lo(asteroid);
        
        if ~isempty(ir_imgs)
            
            % IR image            
            irImg = imread(ir_imgs(j));
            irImg = rgb2gray(irImg);
            irImg(irImg<uint8(20)) = uint8(0);
            irImg = irImg*1000;
            
            % obtain limb
            [trim_u_IR, trim_v_IR,~,~,mid_pt_u_IR,mid_pt_v_IR] = edge_finding_lo(irImg);
        
        end
        %check sign of y comp of SunB
        cam_pos = [0,0,z_list(j)]; %uncomment for regular cases
        sunb = CB(:,:,j)*(cam_pos + sun_pos(j,:))';
        %sunb = cam_pos(j,:) - (CB(:,:,j)'*sun_pos(j,:)')';
        if sunb(2) >= 0
            dir = 1;
        else
            dir = -1;
        end
        
        % optical limb and rays
        [edge_points{j}, edge_rays{j}, new_trim_u,new_trim_v] = edge_to_3d_lo(z_list(j), fov_angle, trim_u, trim_v,sun_pos(j,:),mid_pt_u,mid_pt_v,dir,ext,j);
        
        
        if ~isempty(ir_imgs)
            % ir limb and rays
            [edge_points_IR{j}, edge_rays_IR{j}, new_trim_u_IR,new_trim_v_IR] = edge_to_3d_ir(z_list(j), fov_angle, trim_u_IR, trim_v_IR,sun_pos(j,:),mid_pt_u_IR,mid_pt_v_IR,dir,ext,j);
        end
        if j > 0
            %ast_flip = flip(asteroid,1);
            figure(1)
            imshow(asteroid)
            hold on
            grid on
            scatter(new_trim_u,new_trim_v,'filled','b')
            scatter(mid_pt_u,mid_pt_v,'filled','r')
            if ~isempty(ir_imgs)
                scatter(new_trim_u_IR,new_trim_v_IR,'filled','g')
                legend({'limb','center','IR limb'},'FontSize',24)
            else
                legend({'limb','center'},'FontSize',24)
            end
            %plot lines from sun direction
            %asteroid = imcrop(asteroid);
            xlabel('X (pixels)','FontSize',16)
            ylabel('Y (pixels)','FontSize',16)
            title(string(j),'FontSize',24)
            hold off    
        end
        
        T = CB(:,:,j);% camera to body frame DCM
        
        % convert optical rays and points into body fixed frame
        [edge_points_bc{j}, new_rays_bc] = CFtoBF_orex(T, edge_points{j}, edge_rays{j});

        if ~isempty(ir_imgs)
            % convert IR rays and points into body fixed frame
            [edge_points_bc_IR{j}, new_rays_bc_IR] = CFtoBF_orex(T, edge_points_IR{j}, edge_rays_IR{j});
        end
        
        if isempty(ir_imgs)
            limb_starts{j} = new_rays_bc(:,4:6);
            limb_ends{j}  = new_rays_bc(:,1:3);
        else
            limb_starts{j} = [new_rays_bc(:,4:6);new_rays_bc_IR(:,4:6)];
            limb_ends{j}  = [new_rays_bc(:,1:3);new_rays_bc_IR(:,1:3)];
        end
            
        j = j+1;
    end
end