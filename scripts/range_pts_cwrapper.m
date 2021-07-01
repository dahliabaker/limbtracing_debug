function [range] = range_pts_cwrapper(testPoints, uhat, facets, verts, nhats)
% Function computes the range from each sensor grid point to where it
% intersects the asteroid shape model. 
% Inputs: 
%       tc              current time
%       X               [R_1; V_1; ... R_n; V_n], SF frame 
%       sensor          struct 
%       targets(i)      struct 
%                   *.facets  [n x 3]
%                   *.verts   [m x 3] 
%                   *.nhats   [n x 3] 
%                   *.C_SF_BF
% 
% Outputs:
%       range 
% 
% Ann Dietrich
% 03/27/2016
% based on code from:
% 07/08/2014    original
% 02/06/2014    edit
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        
% Reshape for C 
testPoints = testPoints';
facets = facets';
verts = verts'; 
nhats = nhats';

npts = size(testPoints,2);
nfacets = size(facets,2);
nverts = size(verts,2);

%Reset indices for use in C
facets = facets - 1;


% done in Sensor frame
[range] = range_pts(npts,nfacets,nverts,uhat,testPoints,facets,verts,nhats);

% Reshape outputs
% range = reshape(range,sensor.npix_x,sensor.npix_y);

