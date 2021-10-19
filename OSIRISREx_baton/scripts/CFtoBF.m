 function [coord_BF, rays_BF] = CFtoBF(lat, long, edge_points, edge_rays)

%+90 rotation about first axis
rot1 = [1 0 0; 0 cosd(90) sind(90); 0 -sind(90) cosd(90)];
%-90 rotation about third axis
rot2 = [cosd(270), sind(270), 0; -sind(270), cosd(270), 0; 0, 0, 1];
% -longitude rotation about third axis
rot3 = [cosd(-long) sind(-long) 0; -sind(-long) cosd(-long) 0; 0 0 1];
% -latitude rotation about first axis
rot4 = [1 0 0; 0 cosd(-lat) sind(-lat); 0 -sind(-lat) cosd(-lat)];

BC = rot4*rot3*rot2*rot1;

for i = 1:length(edge_rays)
    rays_BF(i,1:3) = (BC*(edge_rays(i,1:3))')';
    rays_BF(i,4:6) = (BC*(edge_rays(i,4:6))')';

    coord_BF(i,:) = (BC*(edge_points(i,:))')';
end




end