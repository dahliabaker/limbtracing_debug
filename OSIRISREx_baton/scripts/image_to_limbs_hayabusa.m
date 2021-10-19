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

function [limb_starts, limb_ends, edge_points_bc] = image_to_limbs_hayabusa(img_list, z_list, fov_angle,CB,sun_pos,phase,cam_pos,limb,ext)
  
    j = 1;
    
    while j <= length(img_list)

        asteroid = imread(img_list(j));
        %asteroid = rgb2gray(asteroid); %toggle on or off based on input
        asteroid = imadjust(asteroid,[0 0.2]);
        asteroid = imgaussfilt(asteroid,7);
        %imwrite(asteroid,'blur_'+string(j)+'.png')
        

        [trim_u, trim_v,E_u,E_v,mid_pt_u,mid_pt_v] = edge_finding_canny(asteroid,10);

        %check sign of y comp of SunB
        %cam_pos = [z_list(j),0,0];
        sunb = cam_pos(j,:) - (CB(:,:,j)*sun_pos(j,:)')';
        if sunb(2) < 0
            dir = 1;
        else
            dir = -1;
        end
        [edge_points{j}, edge_points_t{j}, edge_rays{j}, edge_rays_t{j}, new_trim_u,new_trim_v, new_term_u,new_term_v] = edge_to_3d_orex(z_list(j), fov_angle, trim_u, trim_v,CB(:,:,j)*sun_pos(j,:)',mid_pt_u,mid_pt_v,dir,phase(j),limb,ext);
        %plot them one over another

        if j ==35
            %ast_flip = flip(asteroid,1);
            figure()
            imshow(asteroid)
            hold on
            scatter(new_trim_u,new_trim_v,'filled','b')
            scatter(new_term_u,new_term_v,'filled','r')
            %plot lines from sun direction
            
            legend({'limb','terminator'},'FontSize',12)
            hold off
            
%             figure()
%             hold on
%             scatter(edge_points{j}(:,1),edge_points{j}(:,2),'b')
%             scatter(edge_points_t{j}(:,1),edge_points_t{j}(:,2),'r')
%             legend({'limb','terminator'},'FontSize',12)
%             hold off
        end
%         saveas(gcf,'hayabusa/output/edge'+string(j)+'.png','png')
  
        %lat = lat_list(j);%degrees
        %long = long_list(j);%degrees
%         h_extra_firstaxis = [1 0 0; 0 cosd(-90) sind(-90); 0 -sind(-90) cosd(-90)]; 
%         h_extra_secondaxis = [cosd(180),0,-sind(180);0,1,0; sind(180),0,cosd(180)];
        T = CB(:,:,j);
   
        sun_v3(j,1:3) = (T*sun_pos(j,:)')./norm(sun_pos(j,:));%-T*cam_pos(j,:)';%in camera frame
        rhat = [0,0,z_list(j)];%camera pointing vector in camera frame
        theta_yz = asind(norm(cross(sun_v3(j,:),rhat))/(norm(sun_v3(j,:))*norm(rhat))); %angle between rhat and sun vector in y-z plane
        %check for which side of camera sun is on
        if sun_v3(j,2) >= 0
            ang = -theta_yz;
        else 
            ang = theta_yz;
        end
        %ang = theta_yz;
        rot1 = [1 0 0; 0 cosd(ang) sind(ang); 0 -sind(ang) cosd(ang)];%ang rotation about first axis
        
%         theta_xz = 180-atan2d(sun_v3(j,1),sun_v3(j,3));
%         rot2 = [cosd(theta_xz) 0 sind(theta_xz); 0 1 0; -sind(theta_xz) 0 cosd(theta_xz)]'; 
        theta_xy = 180-atan2d(sun_v3(j,1),sun_v3(j,2));
        theta_xy = atan2d(sun_v3(j,1),sun_v3(j,2));
        rot2 = [cosd(theta_xy), sind(theta_xy), 0; -sind(theta_xy), cosd(theta_xy), 0; 0, 0, 1]; 
%         
        
         %rot3 = [1 0 0; 0 cosd(90) sind(90); 0 -sind(90) cosd(90)]; 
%         rot4 = [cosd(180),0,-sind(180);0,1,0; sind(180),0,cosd(180)];
        
        rot = rot2*rot1;
        %rot = eye(3,3);
        [edge_points_bc{j}, new_rays_bc] = CFtoBF_orex(T, rot, edge_points{j}, edge_points_t{j}, edge_rays{j},edge_rays_t{j});


%         limb_starts{j} = new_rays_bc(:,1:3);
%         limb_ends{j}  = new_rays_bc(:,4:6);
        limb_starts{j} = new_rays_bc(:,4:6);
        limb_ends{j}  = new_rays_bc(:,1:3);

        
        disp('image number processed: ')
        disp(j)
        
        
        j = j+1;
     end
    

end