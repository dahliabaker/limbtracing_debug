 grav_truth = table2array(readtable('Gravity_Itokawa_truth.csv'));

% Based on your discretization
nlon    = 80;
nlat    = 20;
nL      = 4;

R0      = 0.5; % km
d       = 0.2; % amount of km added for each new radius

lat     = linspace(-pi/2,pi/2,nlat);
lon     = linspace(-pi,pi,nlon);

k = 1;

for i = 1:nL % radii
        
    gravnorm_truth = zeros(nlat,nlon);

    R = R0 + d*(i-1);
    
        for m = 1:nlat % latitudes

            for j = 1:nlon % longitudes

            % put it in matrix form for contour plotting
            gravnorm_truth(m,j) = norm(grav_truth(k,1:3));
            
            k = k + 1;

            end



        end

        subplot(nL/2,nL/2,i)
        hold on
        contourf(rad2deg(lon),rad2deg(lat),gravnorm_truth,50,'LineStyle','none')
        c = colorbar;
        c.Label.String = 'Gravity norm';
        ylabel('Latitude (deg)','Interpreter','Latex')
        xlabel('Longitude (deg)','Interpreter','Latex')
        colormap('jet')
        xticks(-180:90:180)
        yticks(-90:30:90)
        axis equal
        title(['$g$ for $R$ = ' num2str(R) 'km'],'Interpreter','Latex')
        set(gca,'FontSize',15)
        
end