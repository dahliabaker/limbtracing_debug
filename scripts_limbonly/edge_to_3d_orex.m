function [edge_points_woc, edge_points_woc_t, edge_rays, edge_rays_t, new_trim_u,new_trim_v, new_term_u,new_term_v] = edge_to_3d_orex(z, fov_angle, trim_u, trim_v,sun_v,mid_pt_u,mid_pt_v,dir,phase,limb,ext,vis_limb)

edge_points_woc = [];
edge_points_woc_t = [];
edge_rays = [];
edge_rays_t = [];
new_trim_u = [];
new_trim_v = [];

horz_dist = z*tand(fov_angle);
vert_dist = z*tand(fov_angle);

% dist_offset = (offset_v*(vert_dist/1024));
%convert pixel numbers to distances in km from origin
dist_u = (trim_u.*(horz_dist/1024));
dist_v = (trim_v.*(vert_dist/1024));

mean_u = mean(dist_u);
mean_v = mean(dist_v);%+dist_offset;
%try sun_pos in body frame
sun_v = [-sun_v(3)+z,-sun_v(2)];

%how are we defining the sun vector in this frame?

%compute dot product w 2d image sunvector
j=1;
k=1;
for i = 1:length(dist_u)
    vec = [-dist_u(i),-dir*dist_v(i)];%+dist_offset];
    sun = [sun_v(1),sun_v(2)];
    dot_p = dot(vec,sun);

    if dot_p>=0
        new_dist_u(j) = dist_u(i);
        new_dist_v(j) = dist_v(i);

        new_trim_u(j) = trim_u(i)+mid_pt_u;
        new_trim_v(j) = (1*trim_v(i))+mid_pt_v;
        j = j+1;
       
    end
end

half_len = ext;

if j > 1
    max = 100;
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
    
    new_u(len+1,1) = mean_u;
    new_v(len+1,1) = mean_v;
    new_w(len+1,1) = z;

    edge_points = 50.*[new_u, -new_v, new_w];

    for i = 1:len
        edge_points_woc(i,:) = edge_points(i,:) - edge_points(len+1,:);
    end

    %limb_starts
    edge_rays(:,1:3) = edge_points_woc(:,1:3) + [0,0,half_len];
    %limb_ends
    edge_rays(:,4:6) = edge_points_woc(:,1:3) - [0,0,half_len];
    %
end
end