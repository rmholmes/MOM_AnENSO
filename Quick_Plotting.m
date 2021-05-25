
% This script makes forcing files for MOM-AnENSO runs.

clear all;
close all;

addpath(genpath('/short/e14/rmh561/software/matlab-utilities/'));
startup;

% Read from MOM:
base = '/g/data/e14/rmh561/MOM_AnENSO/';
lon = ncread([base 'temp.cat.diff_yr1.nc'],'xt_ocean');
lat = ncread([base 'temp.cat.diff_yr1.nc'],'yt_ocean');
depth = ncread([base 'temp.cat.diff_yr1.nc'],'st_ocean');

[tmp ind1] = min(abs(lon+160));
[tmp ind2] = min(abs(lon+50));
[tmp ind3] = min(abs(lat+78));
[tmp ind4] = min(abs(lat+60));

[tmp dind] = min(abs(depth-200));

time = [ncread([base 'temp.cat.diff_yr1.nc'],'time'); ...
        ncread([base 'temp.cat.diff_yr2.nc'],'time'); ...
        ncread([base 'temp.cat.diff_yr3.nc'],'time')];
tL = length(time);
temp = cat(3,squeeze(ncread([base 'temp.cat.diff_yr1.nc'],'temp',[ind1 ind3 dind 1],[ind2-ind1+1 ind4-ind3+1 1 tL/3])), ...
        squeeze(ncread([base 'temp.cat.diff_yr2.nc'],'temp',[ind1 ind3 dind 1],[ind2-ind1+1 ind4-ind3+1 1 tL/3])), ...
        squeeze(ncread([base 'temp.cat.diff_yr3.nc'],'temp',[ind1 ind3 dind 1],[ind2-ind1+1 ind4-ind3+1 1 tL/3])));

time = time-time(1)+2.5;
[xL, yL] = size(temp(:,:,1));
avgp = 6*3;
temp_3mon = zeros(xL,yL,floor(tL/avgp));
time_3mon = zeros(floor(tL/avgp),1);
for i=1:length(temp_3mon(1,1,:))
    temp_3mon(:,:,i) = mean(temp(:,:,(i-1)*avgp+1:i*avgp),3);
    time_3mon(i) = mean(time((i-1)*avgp+1:i*avgp));
end

[X,Y] = ndgrid(lon(ind1:ind2),lat(ind3:ind4));

amps = zeros(size(time_3mon));
for i=1:length(time_3mon)
    amps(i) = interp1((0.5:1:35.5)/12*365,avgs,time_3mon(i),'linear');
end

months = {'JAS Year 1','OND Year 1','JFM Year 2','AMJ Year 2','JAS Year 2','OND Year 2'};
figure;
for i=3:8
    subplot(2,3,i-2)
    contourf(X,Y,temp_3mon(:,:,i),[-1:0.025:1],'linestyle','none');
    hold on;
    quiver(CX(1:3:end,:)-360,CY(1:3:end,:),U10regCNYF(1:3:end,:)*amps(i)*6.1111,V10regCNYF(1:3:end,:)*amps(i),'-k','AutoScale','off');
    caxis([-0.5 0.5]);
    title([months{i-2} ' 200m T Anomaly ($^\circ$C)']);
    set(gca,'color','k');
    xlabel('Longitude ($^\circ$E)');
    ylabel('Latitude ($^\circ$N)');
    colorbar;
    xlim([-160 -50]);
    ylim([-78 -60]);
end
colormap(redblue);


