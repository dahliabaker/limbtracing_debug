function nhats = ComputeNormals(verts, facets) 
% function to compute the normal vectors for inputed verts and facets
% (normals are normalized)
% nhats = ComputeNormals(verts, facets) 
%   Inputs 
%       verts       [n x 3] 
%       facets      [n x 3]
% 
% Ann Dietrich
% 07/07/2014 


% number of verts
n = size(facets,1); 

% loop through verts and facets 
nhats = NaN(n,3); 
for i = 1:n 
    
    % get indicies from facets 
    facetc = facets(i,:); 
    
    % define needed verts vectors 
    v1 = verts(facetc(1),:); 
    v2 = verts(facetc(2),:);
    v3 = verts(facetc(3),:); 
    
    % compute normal:
    nvec = cross( (v2-v1), (v3-v2) ); 
    % normalize:
    nhats(i,:) = nvec/norm(nvec); 
    
end





