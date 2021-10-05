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

function [limb_starts, limb_ends, edge_points_bc] = image_to_limbs_orex(img_list, z_list, fov_angle,CB,sun_pos,phase,limb,ext,IR,IRdat,density)
  
    j = 1;
    
    while j <= length(img_list)

        asteroid = imread(img_list(j));
        asteroid = rgb2gray(asteroid); %toggle on or off based on input
        %figure for paper
%         asteroid = imadjust(asteroid,[0 1]);
%         asteroid = imgaussfilt(asteroid,5);
%         [x,y]=size(asteroid);
%         X=1:x;
%         Y=1:y;
%         [xx,yy]=meshgrid(Y,X);
%         i=im2double(asteroid);
%         figure;mesh(xx,yy,i);
%         c = colorbar;
%         c.Label.String = 'Pixel Intensity Value';
%         c.Label.FontSize = 12;
%         figure;imshow(i)
%         xlim([0 1024])
%         ylim([0 1024])
%         title('Adjusted and Gaussian Filtered Intensity Map - Itokawa','FontSize',16)
%         ylabel('Vertical Pixel Position','FontSize',12)
%         xlabel('Horizontal Pixel Position','FontSize',12)
%         
%         disp('done')
        %imwrite(asteroid,'blur_'+string(j)+'.png')
        

        
        [trim_u, trim_v,E_u,E_v,mid_pt_u,mid_pt_v] = edge_finding_canny(asteroid, density);
        
        %check sign of y comp of SunB
        cam_pos = [0,0,z_list(j)];
        %sunb = cam_pos - (CB(:,:,j)'*sun_pos(j,:)')';
        sunb = cam_pos - (CB(:,:,j)'*sun_pos(j,:)')';
        if sunb(2) <= 0
            dir = 1;
        else
            dir = -1;
        end
        [edge_points{j}, edge_points_t{j}, edge_rays{j}, edge_rays_t{j}, new_trim_u,new_trim_v, new_term_u,new_term_v] = edge_to_3d_orex(z_list(j), fov_angle, trim_u, trim_v,sunb,mid_pt_u,mid_pt_v,dir,phase,limb,ext);
        %plot them one over another

        if j ==1
            ast_flip = flip(asteroid,1);
            figure()
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
            hold off
            
%             figure()
%             hold on
%             scatter(edge_points{j}(:,1),edge_points{j}(:,2),'b')
%             scatter(edge_points_t{j}(:,1),edge_points_t{j}(:,2),'r')
%             legend({'limb','terminator'},'FontSize',12)
%             hold off
        end
        %saveas(gcf,'output/edge'+string(j)+'.png','png')
  
        %lat = lat_list(j);%degrees
        %long = long_list(j);%degrees
        
        T = CB(:,:,j);
   
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

        disp('image number processed: ')
        disp(j)
        
        
        j = j+1;
    end
    k = 73;
    while k <=(length(IRdat)+72)
        l = k-72;
        [trim_u, trim_v,E_u,E_v,mid_pt_u,mid_pt_v] = edge_finding_IR(IRdat{l,1},IRdat{l,2}, 10);
        
        %check sign of y comp of SunB
        cam_pos = [0,0,z_list(l)];
        %sunb = cam_pos - (CB(:,:,j)'*sun_pos(j,:)')';
        sunb = cam_pos - (CB(:,:,l)'*sun_pos(l,:)')';
        if sunb(2) <= 0
            dir = 1;
        else
            dir = -1;
        end
        [edge_points{k}, edge_points_t{k}, edge_rays{k}, edge_rays_t{k}, new_trim_u,new_trim_v, new_term_u,new_term_v] = edge_to_3d_orex(z_list(l), fov_angle, trim_u, trim_v,sunb,mid_pt_u,mid_pt_v,dir,0,limb,ext);
        %plot them one over another
        
        T = CB(:,:,l);
   
        sun_v3 = sun_pos(l,:)+[0,0,-100];%in camera frame
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
        [edge_points_bc{k}, new_rays_bc] = CFtoBF_orex(T, rot, edge_points{k}, edge_points_t{k}, edge_rays{k},edge_rays_t{k});


        limb_starts{k} = new_rays_bc(:,1:3);
        limb_ends{k}  = new_rays_bc(:,4:6);

        k = k+1;
    end

end