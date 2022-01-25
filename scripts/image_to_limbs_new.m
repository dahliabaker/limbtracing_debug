function [limb_starts, limb_ends, edge_points_bc] = image_to_limbs_new(img_list, z_list, fov_angle,CB,sun_pos,phase,limb,ext,IR,IRdat,density,body)
% initialize cell arrays
limb_starts = cell(length(img_list),1);
limb_ends = cell(length(img_list),1);
edge_points_bc = cell(length(img_list),1);

% add column if IR data exists
if ~isempty(IRdat)
    limb_starts = [limb_starts,limb_starts];
    limb_ends = [limb_ends,limb_ends];
    edge_points_bc = [edge_points_bc,edge_points_bc];
end

% Loop through Both Image sets concurrently
for iImg = 1:length(img_list)
    %% get IR and Optical edges    
    % process optical edge
    im_raw = imread(img_list(iImg));
    asteroid = rgb2gray(im_raw); %toggle on or off based on input
    asteroid(asteroid<uint8(2)) = uint8(0);
    asteroid = asteroid*1000;
    [trim_u_opt, trim_v_opt,E_u_opt,E_v_opt,mid_pt_u,mid_pt_v] = edge_finding_canny(asteroid, density,limb);
    
    % if IR data exists, process that first to get more accurate center
    if ~isempty(IRdat)
        
        [trim_u_IR, trim_v_IR,E_u_IR,E_v_IR,mid_pt_u,mid_pt_v] = edge_finding_IR(IRdat{iImg,2},IRdat{iImg,1}, 5);% change last input to 5 to match optical sampling freq
        
    end
    
    % get sun in body frame
    cam_pos = [0,0,z_list(iImg)];
    sunb = cam_pos - (CB(:,:,iImg)'*sun_pos(iImg,:)')';
    %check sign of y comp of sunb
    if sunb(2) >= 0
        dir = 1;
    else
        dir = -1;
    end
    
    %% process optical edge to get limb and terminator
    [edge_points{iImg}, ...% need not be saved accross iterations
     edge_points_t{iImg}, ...
     edge_rays{iImg}, ...
     edge_rays_t{iImg}, ...
     new_trim_u_opt, ...
     new_trim_v_opt, ...
     new_term_u_opt, ...
     new_term_v_opt] = edge_to_3d_orex(z_list(iImg),...
                                      fov_angle, ...
                                      trim_u_opt, ...
                                      trim_v_opt, ...
                                      sunb, ...
                                      mid_pt_u, ...
                                      mid_pt_v,...
                                      dir,...
                                      phase(iImg), ...
                                      limb, ...
                                      ext);

    % transform optical edge points from camera frame to body frame
    T_c2b = CB(:,:,iImg);
    %in camera frame
    sun_v3 = sun_pos(iImg,:)+[0,0,-100];
    %camera pointing vector in camera frame
    rhat = [0,0,-1];
    %angle between rhat and sun vector in y-z plane
    theta_yz = asind(norm(cross(sun_v3,rhat))/(norm(sun_v3)*norm(rhat))); 
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

    % convert edges from camera fraome to body frame
    [edge_points_bc{iImg,1}, ...
     new_rays_bc] = CFtoBF_orex(T_c2b, ...
     							rot, ...
     							edge_points{iImg}, ...
     							edge_points_t{iImg}, ...
     							edge_rays{iImg}, ...
     							edge_rays_t{iImg});

    limb_starts{iImg,1} = new_rays_bc(:,1:3);
    limb_ends{iImg,1}  = new_rays_bc(:,4:6);

    
    %% process IR edge to get limb ( if IR data exists)
    if ~isempty(IRdat)
        
        [edge_points{iImg}, ...
		 edge_points_t{iImg}, ...
		 edge_rays{iImg}, ...
		 edge_rays_t{iImg}, ...
		 new_trim_u_IR, ...
		 new_trim_v_IR, ...
		 new_term_u_IR, ...
		 new_term_v_IR] = edge_to_3d_orex(z_list(iImg), ...
                                          fov_angle, ...
									      trim_u_IR, ...
									      trim_v_IR, ...
									      sunb, ...
							   		      mid_pt_u, ...
									      mid_pt_v, ...
									      dir, ...
									      0, ...
									      limb, ...
									      ext);

        
    end
        
    % convert edges from camera fraome to body frame
    [edge_points_bc{iImg,2}, ...
     new_rays_bc] = CFtoBF_orex(T_c2b, ...
     							rot, ...
     							edge_points{iImg}, ...
     							edge_points_t{iImg}, ...
     							edge_rays{iImg}, ...
                                edge_rays_t{iImg});
    
    limb_starts{iImg,2} = new_rays_bc(:,1:3);
    limb_ends{iImg,2}  = new_rays_bc(:,4:6);
    fprintf('image % i of % i processed\n',iImg,length(img_list))
    % plotting
    figure(1)
    % imshow(im_raw)
    imshow(asteroid)
    hold on
    % optical data
    scatter(new_trim_u_opt,new_trim_v_opt,'filled','b')% limb
    scatter(new_term_u_opt,new_term_v_opt,'filled','r')% terminator
    
    % ir data
    scatter(new_trim_u_IR,new_trim_v_IR,'filled','g')% limb
    legend({'optical limb','terminator','IR limb'},'FontSize',24)
    xlabel('X (pixels)','FontSize',16)
    ylabel('Y (pixels)','FontSize',16)
    hold off    
    
end
% vectorize cell arrays
limb_starts = limb_starts(:);
limb_ends = limb_ends(:);
edge_points_bc = edge_points_bc(:);