function [coord_BF, rays_BF] = CFtoBF_orex(CB, R, edge_points,edge_points_t, edge_rays,edge_rays_t)
rays_BF = [];
coord_BF = [];
rays_BF_t = [];
coord_BF_t = [];

BC = CB';
%incorporate a translation between the frames
%points first
% edge_points(:,3) = edge_points(:,3) + 100;
% %now rays
% edge_rays(:,3) = edge_rays(:,3) + 100;
% edge_rays(:,6) = edge_rays(:,6) + 100;

for i = 1:length(edge_rays)
    rays_BF(i,1:3) = (BC*(edge_rays(i,1:3))')';
    rays_BF(i,4:6) = (BC*(edge_rays(i,4:6))')';

    coord_BF(i,:) = (BC*(edge_points(i,:))')';
    
end

for i = 1:length(edge_rays_t)
    rays_BF_t(i,1:3) = (BC*R*(edge_rays_t(i,1:3))')';
    rays_BF_t(i,4:6) = (BC*R*(edge_rays_t(i,4:6))')';
    
    coord_BF_t(i,:) = (BC*R*(edge_points_t(i,:))')';
    
end

rays_BF = [rays_BF;rays_BF_t];
coord_BF = [coord_BF; coord_BF_t];

end
