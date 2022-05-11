function [intLineSeg] = patch_patch_intersection(patch1,patch2)

% compute patch normals
n1 = cross(patch1(:,1)-patch1(:,2),...
           patch1(:,3)-patch1(:,2));n1 = normc(n1);
n2 = cross(patch2(:,1)-patch2(:,2),...
           patch2(:,3)-patch2(:,2));n2 = normc(n2);
       
% check if patches are parallel
if norm(cross(n1,n2))==0
    
    intLineSeg = [];
    
else% if not parallel, intersection will exist
    
    % intersection line parameters
    h1 = dot(n1,patch1(:,2));
    h2 = dot(n2,patch2(:,2));
    
    c1 = (h1 - h2*dot(n1,n2))/(1-dot(n1,n2))^2;
    c2 = (h2 - h1*dot(n1,n2))/(1-dot(n1,n2))^2;
    
    
    
end