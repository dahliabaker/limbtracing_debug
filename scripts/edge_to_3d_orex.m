function [edge_points_woc, edge_points_woc_t, edge_rays, edge_rays_t, new_trim_u,new_trim_v, new_term_u,new_term_v] = edge_to_3d_orex(z, fov_angle, trim_u, trim_v,sun_v,mid_pt_u,mid_pt_v,dir,phase,limb,ext)

edge_points_woc = [];
edge_points_woc_t = [];
edge_rays = [];
edge_rays_t = [];
new_trim_u = [];
new_trim_v = [];
new_term_u = [];
new_term_v = [];

%3 polynomial fit
%offset_v = (-3.227e-5*(phase^3)+0.0079*(phase^2)+0.2613*phase - 2.6469);

%4 polynomial fit
offset_v = dir*(3.2966e-7*(phase^4) - 1.3776e-4*(phase^3) + 0.0187*(phase^2) -0.1074*(phase));

%fixing mid_pt_v for flip from origin in top left corner
%mid_pt_v = 512-(512-mid_pt_v);
%Put them in 3D space
% w is going to be the third image dimension

%z = (focal_len*known_width)/app_width;

horz_dist = z*tand(fov_angle);
vert_dist = z*tand(fov_angle);

dist_offset = (offset_v*(vert_dist/1024));
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
    vec = [dist_u(i),-dir*dist_v(i)];%+dist_offset];
    sun = [-sun_v(1),sun_v(2)];
    dot_p = dot(vec,sun);
    
    if (dot_p<=0)
        new_dist_u(j) = dist_u(i);
        new_dist_v(j) = dist_v(i);
        
        new_trim_u(j) = trim_u(i)+mid_pt_u;
        new_trim_v(j) = (1*trim_v(i))+mid_pt_v;
        j = j+1;
    elseif dot_p > 0 && limb == 1 
        term_dist_u(k) = dist_u(i);
        term_dist_v(k) = dist_v(i);
        
        new_term_u(k) = trim_u(i)+mid_pt_u;
        new_term_v(k) = (1*trim_v(i))+mid_pt_v;
        k = k+1;
    end
end

half_len = ext;

if j > 1
    dist_u = new_dist_u';
    dist_v = new_dist_v';
    %Make rays for limb
    len = length(dist_u(:,1));
    dist_w(1:len,1) = z;
    
    dist_u(len+1,1) = mean_u;
    dist_v(len+1,1) = mean_v;
    dist_w(len+1,1) = z;

    edge_points = [dist_u, -dist_v, dist_w];

    for i = 1:len
        edge_points_woc(i,:) = edge_points(i,:) - edge_points(len+1,:);
        mx(i) = edge_points_woc(i,1)/z;
        my(i) = edge_points_woc(i,2)/z;
    end

    %calculate delta x and delta y for slope at half_len distance from pt
    %limb_starts
    %edge_rays(:,1:3) = [(edge_points_woc(:,1) + (mx(:)*half_len)), (edge_points_woc(:,2)+ (my(:)*half_len)),(edge_points_woc(:,3) +half_len)];
    edge_rays(:,1:3) = edge_points_woc(:,1:3) + [0,0,half_len];
    %limb_ends
    %edge_rays(:,4:6) = [(edge_points_woc(:,1) - (mx(:)*half_len)), (edge_points_woc(:,2)- (my(:)*half_len)),(edge_points_woc(:,3) -half_len)];
    edge_rays(:,4:6) = edge_points_woc(:,1:3) - [0,0,half_len];
    %
end

if k > 1
    term_u = term_dist_u';
    term_v = term_dist_v';
    
    len_t = length(term_u(:,1));
    term_w(1:len_t,1) = z;

    term_u(len_t+1,1) = mean_u;
    term_v(len_t+1,1) = mean_v;
    term_w(len_t+1,1) = z;
    
    edge_points_t = [term_u, -term_v, term_w];
    
    for i = 1:len_t
        edge_points_woc_t(i,:) = edge_points_t(i,:) - edge_points_t(len_t+1,:);
        mx(i) = edge_points_woc_t(i,1)/z;
        my(i) = edge_points_woc_t(i,2)/z;
    end
    
    %limb_starts
    %edge_rays_t(:,1:3) = [(edge_points_woc_t(:,1) + (mx(:)*half_len)), (edge_points_woc_t(:,2) + (my(:)*half_len)),(edge_points_woc_t(:,3) +half_len)];
    edge_rays_t(:,1:3) = edge_points_woc_t(:,1:3) + [0,0,half_len];
    %edge_rays_t(:,1:3) = edge_points_woc_t(:,1:3) + [half_len,0,0];
    %limb_ends
    %edge_rays_t(:,4:6) = [(edge_points_woc_t(:,1) - (mx(:)*half_len)), (edge_points_woc_t(:,2) - (my(:)*half_len)),(edge_points_woc_t(:,3) -half_len)];
    edge_rays_t(:,4:6) = edge_points_woc_t(:,1:3) - [0,0,half_len];
    %edge_rays_t(:,4:6) = edge_points_woc_t(:,1:3) - [half_len,0,0];
end

end