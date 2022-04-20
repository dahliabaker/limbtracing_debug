function [edge_points_woc, edge_points_woc_t, edge_rays, edge_rays_t, new_trim_u,new_trim_v,new_term_u,new_term_v] = edge_to_3d_lo(z, fov_angle, trim_u, trim_v,sun_v,mid_pt_u,mid_pt_v,dir,ext,img_num)

edge_points_woc = [];
edge_points_woc_t = [];
edge_rays = [];
edge_rays_t = [];
new_trim_u = [];
new_trim_v = [];
new_term_u = [];
new_term_v = [];
horz_dist = z*tand(fov_angle);
vert_dist = z*tand(fov_angle);

% dist_offset = (offset_v*(vert_dist/1024));
%convert pixel numbers to distances in km from origin
dist_u = (trim_u.*(horz_dist/1024));
dist_v = (trim_v.*(vert_dist/1024));

mean_u = mean(dist_u);
mean_v = mean(dist_v);%+dist_offset;
%try sun_pos in body frame
%sun_v = [-sun_v(3)+z,-sun_v(2)];

%how are we defining the sun vector in this frame?

%compute dot product w 2d image sunvector
j=1;
k=1;
for i = 1:length(dist_u)
    vec = [dist_u(i),-dist_v(i)];%+dist_offset];
    sun = [sun_v(1),sun_v(2)];
    dot_p = dot(vec,sun);
    dot_param = 0;
%     if (img_num <15 && img_num>3) || (img_num > 40 && img_num < 50)
%         dot_param = -15;
%     else
%         dot_param = 0;
%     end
    
    if dot_p<=dot_param
        new_dist_u(j) = dist_u(i);
        new_dist_v(j) = dist_v(i);

        new_trim_u(j) = trim_u(i)+mid_pt_u;
        new_trim_v(j) = (1*trim_v(i))+mid_pt_v;
        j = j+1;
       
    else
        term_dist_u(k) = dist_u(i);
        term_dist_v(k) = dist_v(i);

        new_term_u(k) = trim_u(i)+mid_pt_u;
        new_term_v(k) = (1*trim_v(i))+mid_pt_v;
        k = k+1;
    end
end

half_len = ext;

if j > 1
    max = 25;
    %make sample indices
    samp = 1:(length(new_dist_u)/max):length(new_dist_u);
    samp = round(samp);
    for i = 1:max
        new_u(i) = new_dist_u(samp(i));
        new_v(i) = new_dist_v(samp(i));
    end
    
    new_u = new_u';
    new_v = new_v';
    %Make rays for limb
    len = length(new_u(:,1));
    new_w(1:len,1) = z;
    
    new_u(len+1,1) = 0;
    new_v(len+1,1) = 0;
    new_w(len+1,1) = z;

    edge_points = [new_u, -new_v, new_w];

    for i = 1:len
        edge_points_woc(i,:) = edge_points(i,:) - edge_points(len+1,:);
    end

    %limb_starts
    edge_rays(:,1:3) = edge_points_woc(:,1:3) + [0,0,half_len];
    %limb_ends
    edge_rays(:,4:6) = edge_points_woc(:,1:3) - [0,0,half_len];
    %
end
if k > 1
    max = 25;
    %make sample indices
    samp = 1:(length(term_dist_u)/max):length(term_dist_u);
    samp = round(samp);
    for i = 1:max
        term_u_2(i) = term_dist_u(samp(i));
        term_v_2(i) = term_dist_v(samp(i));
    end
    term_u_2 = term_u_2';
    term_v_2 = term_v_2';
    
    len_t = length(term_u_2(:,1));
    term_w_2(1:len_t,1) = z;
    term_u_2(len_t+1,1) = mean_u;
    term_v_2(len_t+1,1) = mean_v;
    term_w_2(len_t+1,1) = z;
    
    edge_points_t = [term_u_2, -term_v_2, term_w_2];
    
    for i = 1:len_t
        edge_points_woc_t(i,:) = edge_points_t(i,:) - edge_points_t(len_t+1,:);
    end
    
    %limb_starts
    edge_rays_t(:,1:3) = edge_points_woc_t(:,1:3) + [0,0,half_len];
    %limb_ends
    edge_rays_t(:,4:6) = edge_points_woc_t(:,1:3) - [0,0,0];
end
end