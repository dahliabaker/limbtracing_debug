function [limbR, AzTest, limbPos, BS, noiseApplied] = create_limbs_equator(shapeFile, long_sensor, numAz, num_bisect, noiseVals)
%
% Inputs:
%
%       noiseVals = 3 numbers, first is sigma for limb radius, second for
%       azimuth, third for longitude: if 0, no noise is applied
%
% Outputs:
%
%       noiseApplied - sampled values applied, cell array, vector for
%       range noises, vector for az noises, one value for longitude noise
%       if no noise was applied, that cell is empty
%
%
%
% Jay mods for RPI limb finding
% Get a list of points to test in Az/R coordinates
% Get range of these points
% bisect last hit and first miss at that Az a set number of times
% return the R value for all Az values requested 


%% Create initial list of points to test
% ======================= 
% Asteroid Shape Model
% ======================= 

% In body-fixed frame
[facets, verts] = read_obj(shapeFile);
nhats = ComputeNormals(verts,facets);


% Assumes body is centered in the sensor frame
R = sqrt(diag(verts*verts'));
Rmax = max(R);
Rmin = min(R);
AzTest = linspace(0,360,numAz+1).*pi./180.0;
AzTest(end) = [];

noiseApplied = cell(1,3);

% add azimuth noise if requested
if noiseVals(2) ~= 0
    aznoise = noiseVals(2).*randn(1,numAz).*pi./180.0;
    noiseApplied{2} = aznoise;
    AzTest = AzTest + aznoise;
end

% List of points to test in Radius, Az space
ptList = zeros(3.*length(AzTest),2);
for ii = 1:length(AzTest)
    ptList(3*ii-2,:) = [.99*Rmin AzTest(ii)];
    ptList(3*ii-1,:) = [(.99*Rmin + 1.01*Rmax)/2.0 AzTest(ii)];
    ptList(3*ii,:) = [1.01*Rmax AzTest(ii)];
end


%% Set up sensor location/frame
% first case will be viewing from equatorial plane
% uhat direction is always in the center radial direction; not a lens model

r_sensor = 2*Rmax;
if noiseVals(3) ~= 0
    longNoise = noiseVals(3).*randn(1).*pi./180.0;
    noiseApplied{3} = longNoise;
    uhat = -[cos(long_sensor+longNoise); sin(long_sensor+longNoise); 0];
else
    uhat = -[cos(long_sensor); sin(long_sensor); 0];
end
sensor_loc = -r_sensor.*uhat;
up_hat = [0;0;1];
right_hat = cross(uhat,up_hat);

% DCM sensor frame to asteroid body frame
BS = [right_hat up_hat uhat];


%% Get initial ray tracing points
testPoints = zeros(length(ptList),3);
for ii = 1:length(ptList)
    testPoints(ii,:) = sensor_loc + ptList(ii,1).*BS*[cos(ptList(ii,2)); sin(ptList(ii,2)); 0];
end

%% Take initial range measurements
[range] = range_pts_cwrapper(testPoints, uhat, facets, verts, nhats);
indices_bisect = zeros(length(AzTest),1); % highest range value that hits target
for jj = 1:length(AzTest)
    if ~isnan(range(3*(jj-1)+2))
        indices_bisect(jj) = 2;
    else
        indices_bisect(jj) = 1;
    end
end

% Get first bisected test point
testPoints_bisect = zeros(length(AzTest),3);
for jj = 1:length(AzTest)
        
    testPoints_bisect(jj,:) = (testPoints(3*(jj-1) + indices_bisect(jj),:) + testPoints(3*(jj-1) + 1 + indices_bisect(jj),:))./2;
    
end


%% Iterate through fixed number of bisections to refine limb locations
for ii = 2:num_bisect
    
    [range] = range_pts_cwrapper(testPoints_bisect, uhat, facets, verts, nhats);
    for jj = 1:length(AzTest)
        
        if isnan(range(jj))
            
            testPoints(3*(jj-1) + 1 + indices_bisect(jj),:) = testPoints_bisect(jj,:);
            testPoints_bisect(jj,:) = (testPoints(3*(jj-1) + indices_bisect(jj),:) + testPoints(3*(jj-1) + 1 + indices_bisect(jj),:))./2;
        
        else
            
            testPoints(3*(jj-1) + indices_bisect(jj),:) = testPoints_bisect(jj,:);
            testPoints_bisect(jj,:) = (testPoints(3*(jj-1) + indices_bisect(jj),:) + testPoints(3*(jj-1) + 1 + indices_bisect(jj),:))./2;
            
        end
    
    end

end

%% Save last hitting points
limbPos = zeros(length(AzTest), 3);
limbR = zeros(length(AzTest), 1);
[range] = range_pts_cwrapper(testPoints_bisect, uhat, facets, verts, nhats);

% add radius noise if requested
if noiseVals(1) ~= 0
    rnoise = noiseVals(1).*randn(numAz,1);
    noiseApplied{1} = rnoise;
end

for jj = 1:length(AzTest)
    
    if isnan(range(jj))
        
        limbPos(jj,:) = testPoints(3*(jj-1) + indices_bisect(jj),:);
        
    else
        
        limbPos(jj,:) =  testPoints_bisect(jj,:);
        
    end
    
    if noiseVals(1) ~= 0
        tempPos = limbPos(jj,:)' - sensor_loc;
        tempPosMag = norm(tempPos);
        tempPos = tempPos + rnoise(jj).*(tempPos./tempPosMag);
        limbPos(jj,:) = (tempPos + sensor_loc)';
    end
    
    limbR(jj) = norm(limbPos(jj,:)' - sensor_loc);
    
end

end