function [m, x0] = two_plane_intersection(n1, n2, d1, d2)
%
% Compute the line that describes the intersection of two planes
% 
% Inputs:
%   n1 and n2 - normal vectors of the two planes
%
% Outputs:
%   m - direction vector of the line
%   b - a point on the line
%
% Ref: http://geomalgorithms.com/a05-_intersect-1.html

% Direction is cross product of two normal vectors

m = cross(n1,n2);
mmag = norm(m);
if mmag < 1e-10
    x0 = NaN.*ones(3,1);
    return
end
m = m./mmag;

% Find a point on the line
[mMax, mInd] = max(abs(m));

A = [n1'; n2'];
A(:,mInd) = [];
d = -1.*[d1; d2];

b = A\d;

x0 = zeros(3,1);
if mInd == 1
    x0(2) = b(1);
    x0(3) = b(2);
elseif mInd == 2
    x0(1) = b(1);
    x0(3) = b(2);
else
    x0(1:2) = b;
end

end