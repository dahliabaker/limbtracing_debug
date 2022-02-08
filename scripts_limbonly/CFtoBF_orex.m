function [coord_BF, rays_BF] = CFtoBF_orex(CB, edge_points, edge_rays)
rays_BF = [];
coord_BF = [];

BC = CB';


%fudge factor to eliminate silhouette overlap
%gonna push terminator by like 0.001 degrees


for i = 1:length(edge_rays)
    rays_BF(i,1:3) = (BC*(edge_rays(i,1:3))')';
    rays_BF(i,4:6) = (BC*(edge_rays(i,4:6))')';

    coord_BF(i,:) = (BC*(edge_points(i,:))')';
    
end

rays_BF = [rays_BF];
coord_BF = [coord_BF];

end
