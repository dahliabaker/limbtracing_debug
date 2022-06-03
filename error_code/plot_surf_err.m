%deal with ply data


%open ply file data
load('itokawa_15_lo.mat');
data = itokawamodello; 

pt_err = zeros(length(data),1);

%calculate lat and long for each point
for i = 1:length(data)
    x = data(i,1);
    y = data(i,2);
    z = data(i,3);
    pt_err(i) = data(i,4)*1000;
    
    lat(i) = atan2d(z,abs(y));
    long(i) = atan2d(y,x);
end

%bin the spherical points somehow to make an even distribution

lat_bin = linspace(-90,90,180);
lon_bin = linspace(-180,180,180);

%sort the latitudes

[xq,yq,vq] = griddata(long,lat,pt_err',lon_bin,lat_bin');


abs_mean = mean(abs(pt_err));


%%

figure()
contourf(xq,yq,vq,'LineStyle','none')
c = colorbar;
c.Label.String = 'Surface Error (m)';
xlim([-180 180])
ylim([-90 90])
ylabel('Latitude (deg)','Interpreter','Latex')
xlabel('Longitude (deg)','Interpreter','Latex')
colormap('jet')
xticks(-180:30:180)
yticks(-90:30:90)
axis equal
title(['Bennu Model Error - $15^{\circ}$ Terminator-Only'],'Interpreter','Latex')
set(gca,'FontSize',15)
        


%%
%regular histograms

figure()
histogram(pt_err)
grid on
ylabel('Frequency','Interpreter','Latex')
xlabel('Vertex Error (m)','Interpreter','Latex')
title(['Vertex Error Distribution - Itokawa Limb-Only, $15^{\circ}$ Phase'],'Interpreter','Latex')
set(gca,'FontSize',15)


%%

%collapse data along longitude to get a mean for a lat
for i = 2:179
    lat_avg_2(i) = mean(vq(2:179,i),'omitnan');
end



%%
load('i_lo.mat')
lat_bin = linspace(-90,90,180);
%fancy plots

figure()
hold on
plot(lat_bin(1:179),lat_avg_1,'-x')
plot(lat_bin(1:179),lat_avg_2,'-x')
%plot(lat_bin(1:179),lat_avg_3,'-x')
grid on
xlim([-90 90])
xlabel('Latitude (degrees)','Interpreter','Latex')
ylabel('Mean Error (m)','Interpreter','Latex')
legend('$0^{\circ}$','$15^{\circ}$','Interpreter','Latex','FontSize',16)
title(['Error by Latitude - Itokawa Limb-Only Models'],'Interpreter','Latex')
set(gca,'FontSize',18)


