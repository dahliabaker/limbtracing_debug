function [shapePnts, t_intersect, line_uhat, line_x0, nhat1, nhat2] = find_shape_pts_from_limb_segments(limb1, limb2, numPnts)
%
% Test to see if two limb segments' planes intersect, and if they do
% return points along those segments
%
% Each limb section plane is defined by 4 corners in the inputs here, 
% passed in as a 3 x 4 matrix, the first two points or the
% first two columns are the sensor locations used to generate in the order
% they were generated, column 3 is the end of the column 1 ray, and 4 with 2
%
% numPnts - maximum number of points on the body that will be returned
%
% Assumes a rectangular limb plane segment


% assign empty outputs in case of early return
shapePnts = [];
t_intersect = [];


[nhat1, d1, PB1, s1mag1, s2mag1] = get_plane(limb1);
[nhat2, d2, PB2, s1mag2, s2mag2] = get_plane(limb2);


% Find intersection line of two planes, if it exists
[line_uhat, line_x0] = two_plane_intersection(nhat1, nhat2, d1, d2);

if isnan(line_x0(1))
    return
end

% Check to see if line falls within plane segments
% First plane computes set of points if within
t1 = check_plane_segment(PB1, s1mag1, s2mag1, limb1(:,1), line_uhat, line_x0);
if isempty(t1)
    return
end

% Second plane just tests these points
% If any left, will be returned, otherwise an empty matrix will be returned
% shapePnts = test_points_in_plane_segment(PB2, s1mag2, s2mag2, limb2(:,1), shapePnts);
t2 = check_plane_segment(PB2, s1mag2, s2mag2, limb2(:,1), line_uhat, line_x0);
if isempty(t2)
    return
end

% Keep the intersection of the arcs that intersects both patches
t_intersect = [max([t1(1) t2(1)]) min([t1(2) t2(2)])];

if t_intersect(1) <= t_intersect(2)
    shapePnts = linspace(t_intersect(1),t_intersect(2),numPnts).*line_uhat + line_x0;
else
    shapePnts = [];
end

end

% sub-functions
function [nhat, d, PB, s1mag, s2mag] = get_plane(limbMat)
% based on my assumptions, this will give an outward normal I think...
% can test this by making sure d is positive
s1 = limbMat(:,3) - limbMat(:,1);
s1mag = norm(s1);
s1 = s1./s1mag;
s2 = limbMat(:,2) - limbMat(:,1);
s2mag = norm(s2);
s2 = s2./s2mag;
nhat = cross(s1,s2);
nhat = nhat./norm(nhat);

% Body-to-plane rotation matrix
PB = [s1 s2 nhat]';

d = -dot(nhat, limbMat(:,1));

end

function [tCross] = check_plane_segment(PB, s1mag, s2mag, Porigin, uhat, x0)
% See if the intersecting line crosses the segment of the plane we care about

% Convert line equations to plane frame
uhatP = PB*uhat;
x0P = PB*(x0 - Porigin);

% uhatP(3) and x0P(3) should be zero!

% check for special cases of vertical/horizontal slopes
if (abs(uhatP(1)) > 1e-12) && (abs(uhatP(2)) > 1e-12)
    
    % get regular form of line equation in 2D
    
    m = uhatP(2)/uhatP(1);
    b = x0P(2) - m*x0P(1);
    
    xAxis_intercept = b;                % y value when x = 0
    yAxis_intercept = -b/m;             % x value when y = 0
    patchX_intercept = m*s1mag + b;     % y value when x = s1mag
    patchY_intercept = (s2mag-b)/m;     % x value when y = s2mag
    
    tCross = [];
    if (xAxis_intercept >= 0 && xAxis_intercept <= s2mag)
        tCross = -x0P(1)/uhatP(1);
    end
    if (yAxis_intercept >= 0 && yAxis_intercept <= s1mag)
        tCross = [tCross; -x0P(2)/uhatP(2)];
    end
    if (patchX_intercept >= 0 && patchX_intercept <= s2mag)
        tCross = [tCross; (s1mag-x0P(1))/uhatP(1)];
    end
    if (patchY_intercept >= 0 && patchY_intercept <= s1mag)
        tCross = [tCross; (s2mag-x0P(2))/uhatP(2)];
    end
    
    
%     if length(tCross) == 2 % intersects!
%         intersectPnts = linspace(tCross(1),tCross(2),numPnts).*uhat + x0;
%     else
    if length(tCross) < 2 % no intersects!
        tCross = [];
    end
    
elseif abs(uhatP(1) <= 1e-12) % "vertical" line
    
    if x0P(1) >= 0 && x0P(1) <= s1mag
        tCross = [-x0P(2)/uhatP(2); (s2mag-x0P(2))/uhatP(2)];
    else
        tCross = [];
    end
    
elseif abs(uhatP(2) <= 1e-12) % "horizontal" line
    
    if x0P(2) >= 0 && x0P(2) <= s2mag
        tCross = [-x0P(1)/uhatP(1); (s1mag-x0P(1))/uhatP(1)];
    else
        tCross = [];
    end
    
end

tCross = sort(tCross); 

end

function [intersectPnts] = test_points_in_plane_segment(PB, s1mag, s2mag, Porigin, testPnts)

testPntsP = PB*(testPnts - Porigin);
pntsMask = true(1, size(testPntsP,2));
for ii = 1:size(testPntsP,2)
    if (testPntsP(1,ii) < 0 || testPntsP(1,ii) > s1mag || testPntsP(2,ii) < 0 || testPntsP(2,ii) > s2mag)
        pntsMask(ii) = false;
    end
end

intersectPnts = testPnts(:,pntsMask);

end