function [edge_points_woc, edge_points_woc_t, edge_rays, edge_rays_t, new_trim_u,new_trim_v,new_term_u,new_term_v] = edge_to_3d_term(z, fov_angle, trim_u, trim_v,trim_u_ir,trim_v_ir,sun_v,mid_pt_u,mid_pt_v,dir,ext,img_num)

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
%%
% ax = subplot(121);cla

% ir_im = rgb2gray(imread(ir_imgs(img_num)));
% ir_im(ir_im<uint8(2)) = 0;
% ir_edge = edge(ir_im*1000,'canny',.6,10);
% imshow(ir_edge), hold on
% [trim_v_ir,trim_u_ir] = find(ir_edge);
% scatter(trim_u_ir,trim_v_ir,'filled','r'), hold on
% scatter(trim_u+mid_pt_u,trim_v+mid_pt_v,'filled','g'), hold on,axis equal
pxTol = 10;
% for ii = 1:length(trim_u)
%    
%     pxInd = [trim_u(ii),trim_v(ii)]+[mid_pt_u,mid_pt_v];
%     
%     % compare against ir pixels
%     pxDiff = vecnorm(pxInd-[trim_v_ir,trim_u_ir],2,2);
%     if min(pxDiff) < pxTol% if close to ir pixel, then ID as limb
%         scatter(pxInd(1),pxInd(2),'filled','b')
%     else
%         scatter(pxInd(1),pxInd(2),'filled','w')
%     end
%     drawnow
% end

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
%     
%     if dot_p<=dot_param
%         new_dist_u(j) = dist_u(i);
%         new_dist_v(j) = dist_v(i);
% 
%         new_trim_u(j) = trim_u(i)+mid_pt_u;
%         new_trim_v(j) = (1*trim_v(i))+mid_pt_v;
%         j = j+1;
%        
%     else
%         term_dist_u(k) = dist_u(i);
%         term_dist_v(k) = dist_v(i);
% 
%         new_term_u(k) = trim_u(i)+mid_pt_u;
%         new_term_v(k) = (1*trim_v(i))+mid_pt_v;
%         k = k+1;
%     end
    pxInd = [trim_u(i);trim_v(i)]+[mid_pt_u;mid_pt_v]; 
    ir_px = [trim_u_ir;trim_v_ir]+[mid_pt_u;mid_pt_v]; 
    
    % compare against ir pixels
    pxDiff = vecnorm(pxInd-ir_px);
    if min(pxDiff) < pxTol% if close to ir pixel, then ID as limb
        new_dist_u(j) = dist_u(i);
        new_dist_v(j) = dist_v(i);
        
        if trim_v(i)<0
            new_trim_u(j) = trim_u(i)+mid_pt_u;
            new_trim_v(j) = (1*trim_v(i))+mid_pt_v;
            j = j+1;
        else
            new_term_u(k) = trim_u(i)+mid_pt_u;
            new_term_v(k) = (1*trim_v(i))+mid_pt_v;
            k = k+1;
        end
%         scatter(trim_u(i)+mid_pt_u,trim_v(i)+mid_pt_v,'filled','k')
    else% else, ID as terminator
        term_dist_u(k) = dist_u(i);
        term_dist_v(k) = dist_v(i);

        new_term_u(k) = trim_u(i)+mid_pt_u;
        new_term_v(k) = (1*trim_v(i))+mid_pt_v;
        k = k+1;
%         scatter(trim_u(i)+mid_pt_u,trim_v(i)+mid_pt_v,'filled','k')
    end
    
end

half_len = ext;

if j > 1
    maxPts = 25;
    %make sample indices
    samp = 1:(length(new_dist_u)/maxPts):length(new_dist_u);
    samp = round(samp);
    for i = 1:maxPts
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
    maxPts = 25;
    %make sample indices
    samp = 1:(length(term_dist_u)/maxPts):length(term_dist_u);
    samp = round(samp);
    for i = 1:maxPts
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