%%
% Update
% 2/3/2022
% Dahlia Baker





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

function [limb_starts, ...
          limb_ends, ...
          edge_points_bc_l, ...
          termRays_starts, ...
          termRays_ends, ...
          edge_points_bc_t] = image_to_limbs_term(img_list, ...
                                                  z_list, ...
                                                  fov_angle, ...
                                                  CB,sun_pos,ext)
  
    % initialize arrays
    edge_points_bc_l = cell(1,length(img_list));
    edge_points_bc_t = cell(1,length(img_list));
    limb_starts      = cell(1,length(img_list));
    limb_ends        = cell(1,length(img_list));
    termRays_starts  = cell(1,length(img_list));
    termRays_ends    = cell(1,length(img_list));
    
    figure(1),clf
    ax = subplot(122); hold on
    % load ir images and reference shape model
    if contains(img_list(1),'itokawa')
        load('ir_list_i.mat','ir_imgs')
        load('./../PhD_work/surfaceGeneration/Data/itokawa200kData.mat','obj')
    else
        load('ir_list_b.mat','ir_imgs')
        load('./../PhD_work/surfaceGeneration/Data/bennu200kData.mat','obj')
    end
    p = patch('Faces',obj.f.v,'Vertices',obj.v);axis(ax,'equal')
    p.FaceColor = 'w';
    p.EdgeColor = [.01,.01,.01];
    j = 1;
    while j <= length(img_list)

        im_raw = imread(img_list(j));
        asteroid = rgb2gray(im_raw); 
%         asteroid(asteroid<uint8(20)) = uint8(0);
        asteroid = asteroid*1000;
          
        [trim_u, trim_v,~,~,mid_pt_u,mid_pt_v] = edge_finding_lo(asteroid);
        
        % do edge finding on ir image
        ir_image_raw = imread(ir_imgs(j));
        ir_image = rgb2gray(ir_image_raw);
%         ir_image(ir_image<uint8(2)) = 0;
        ir_image = ir_image*1000;
        [trim_u_ir, trim_v_ir,~,~,~,~] = edge_finding_lo(ir_image);
        
        %check sign of y comp of SunB
        cam_pos = [0,0,z_list(j)]; %uncomment for regular cases
        sunb = CB(:,:,j)*(cam_pos + sun_pos(j,:))';
        %sunb = cam_pos(j,:) - (CB(:,:,j)'*sun_pos(j,:)')';
        if sunb(2) >= 0
            dir = 1;
        else
            dir = -1;
        end
        
        [edge_points, ...
         edge_points_t, ...
         edge_rays, ...
         edge_rays_t, ...
         new_trim_u, ...
         new_trim_v, ...
         new_term_u, ...
         new_term_v] = edge_to_3d_term(z_list(j), ...
                                       fov_angle, ...
                                       trim_u, ...
                                       trim_v, ...
                                       trim_u_ir, ...
                                       trim_v_ir, ...
                                       sun_pos(j,:), ...
                                       mid_pt_u, ...
                                       mid_pt_v, ...
                                       dir,ext,j);
        %plot them one over another
        if j > 0
            %ast_flip = flip(asteroid,1);
%             figure(1)
            subplot(121)
%             imshow(asteroid)
            imshow(im_raw)
            hold on
            grid on
            scatter(new_trim_u,new_trim_v,'filled','b')
            scatter(new_term_u,new_term_v,'filled','g')
            scatter(mid_pt_v,mid_pt_v,'filled','r')
            if exist('trim_u_ir','var')
%                 scatter(trim_u_ir+mid_pt_u,trim_v_ir+mid_pt_v,'filled','w')
%                 legend({'limb','terminator','center','ir limb'},'FontSize',24)
            else
                legend({'limb','terminator','center'},'FontSize',24)
            end
            xlabel('X (pixels)','FontSize',16)
            ylabel('Y (pixels)','FontSize',16)
            title(string(j),'FontSize',24)
            hold off    
        end
        
        % body to camera frame rotation 
        T_b2c = CB(:,:,j);
        [edge_points_bc_l{j}, ...
         new_rays_bc_l, ....
         edge_points_bc_t{j}, ...
         rays_BF_t] = CFtoBF_term(T_b2c, ...
                                  eye(3,3), ...
                                  edge_points, ...
                                  edge_points_t, ...
                                  edge_rays, ...
                                  edge_rays_t);
        % limb ray start and ends
        limb_starts{j} = new_rays_bc_l(:,1:3);
        limb_ends{j}  = new_rays_bc_l(:,4:6);
       
        % terminator ray start and ends
        termRays_starts{j} = rays_BF_t(:,1:3);
        termRays_ends{j} = rays_BF_t(:,4:6);
        for ii = 1:length(rays_BF_t)
            termRay =[rays_BF_t(ii,1:3);
                      rays_BF_t(ii,4:6)];
            plot3(ax,termRay(:,1),termRay(:,2),termRay(:,3),'g*-.','LineWidth',2)
        end
%         view(ax,90-((j-1)*5),0)
        view(ax,90,90)
%         if 90-((j-1)*5) >0
%             camup(ax,[0 1 0])
%         else
%             camup(ax,[0 -1 0])
%         end
        drawnow, delete(ax.Children(1:end-1))
        % increment counter
        j = j+1;
        
    end
end