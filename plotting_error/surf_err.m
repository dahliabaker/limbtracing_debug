%error surface map plotting


%4/14/22
%Dahlia Baker

load('model_error.mat');
addpath("geodetic299/geodetic/")
%%
for i = 1:length(bennumodellohd1(:,1))
    
    [lat(i),lon(i),h(i)]=xyz2ell3(bennumodellohd1(i,1),bennumodellohd1(i,2),bennumodellohd1(i,3),.278,.278,1);

    
end

[LAT,LONG] = meshgrid(lat,long);