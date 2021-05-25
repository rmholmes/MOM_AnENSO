
% This script plots SLP and wind vector anomalies related to ENSO
% in the Antarctic Region.

clear all;
close all;

addpath(genpath('/short/e14/rmh561/software/matlab-utilities/'));
startup;

load coastlines
coastlon(coastlon<0) = coastlon(coastlon<0)+360;
inds = [];
for i=2:length(coastlon)
    if (abs(coastlon(i)-coastlon(i-1))>50)
        inds = [inds i];
    end
end
for i=1:length(inds)
    coastlon((inds(i)+1):end+1) = coastlon(inds(i):end);
    coastlat((inds(i)+1):end+1) = coastlat(inds(i):end);
    coastlon(inds(i)+i) = NaN;
    coastlat(inds(i)+i) = NaN;
end
coastlon(372) = NaN;
coastlat(372) = NaN;

type = 1; % 0 = ERA, 1 = JRA
if (type)
    base = 'JRAdata/';
    label = 'JRA55-do';
else
    base = 'ERAdata/';
    label = 'ERA Interim';
end

% Nino 3.4:
DATA = load('index_data/NinoIndices.txt');
n34 = DATA(:,10);
n34yr = DATA(:,1);
n34mn = DATA(:,2);
n34yrd = n34yr + n34mn/12;

% TPI:
IPODATA = load('index_data/tpi.timeseries.ersstv5.filt.dataraw');
iponyrs = length(IPODATA(:,1));
ipo = reshape(IPODATA(:,2:end)',[iponyrs*12 1]);
ipo(abs(ipo)>2) = NaN;
ipoyr = reshape(repmat(IPODATA(:,1),[1 12])',[iponyrs*12 1]);
ipomn = repmat([1:12]',[iponyrs 1]);
ipoyrd = ipoyr + ipomn/12;

% $$$ % ERSST:
% $$$ DATA = load('index_data/nino34.data_ERSSTraw');
% $$$ % HadISST:
% $$$ DATA = load('index_data/nino34.long.data_HadISSTraw');
% $$$ n34nyrs = length(DATA(:,1));
% $$$ n34 = reshape(DATA(:,2:end)',[n34nyrs*12 1]);
% $$$ n34(abs(n34)>50) = NaN;
% $$$ n34yr = reshape(repmat(DATA(:,1),[1 12])',[n34nyrs*12 1]);
% $$$ n34mn = repmat([1:12]',[n34nyrs 1]);
% $$$ n34yrd = n34yr + n34mn/12;
% $$$ n34cli = zeros(12,1);
% $$$ for mi = 1:12
% $$$     n34cli(mi) = nanmean(n34(n34mn == mi & (n34yr >= 1981 & n34yr <=2010)));
% $$$     n34(n34mn == mi) = n34(n34mn == mi) - n34cli(mi);
% $$$ end

Uname = [base 'U10_anom.nc'];
Vname = [base 'V10_anom.nc'];
if (type)
    Pname = [base 'slp_anom.nc'];
else
    Pname = [base 'MSL_anom.nc'];
end
if (type)
    U10 = ncread(Uname,'uas_10m');
    V10 = ncread(Vname,'vas_10m');
    MSL = ncread(Pname,'psl')/100;
    lat = ncread(Uname,'latitude');
    lon = ncread(Uname,'longitude');
    time = ncread(Uname,'time');
    dnum = datenum([1900 1 1 0 0 0])+time;
else
    U10 = ncread(Uname,'10U_GDS0_SFC');
    V10 = ncread(Vname,'10V_GDS0_SFC');
    MSL = ncread(Pname,'MSL_GDS0_SFC')/100;
    lat = ncread(Uname,'g0_lat_1');
    lon = ncread(Uname,'g0_lon_2');
    time = ncread(Uname,'initial_time0_hours');
    dnum = datenum([1800 1 1 0 0 0])+time/24;
end
[xL,yL,tL] = size(U10);

[X,Y] = ndgrid(lon,lat);
dvec = datevec(dnum);
yr = dvec(:,1);
mn = dvec(:,2);
yrd = yr + mn/12;

%% Make an ASL time series:
[tmp ln1] = min(abs(lon-200));
[tmp ln2] = min(abs(lon-300));
[tmp lt1] = min(abs(lat+80));
[tmp lt2] = min(abs(lat+45));
for ti=1:tL
    ASL(ti) = min(min(MSL(ln1:ln2,lt1:lt2,ti)));
end


%%%%%% Calculate N34 regression:
minyr = max([min(yr) min(n34yr)]);
maxyr = min([max(yr) max(n34yr)]);
minyr = 1982;

n34i = find(n34yr==minyr,1,'first');
n34f = find(n34yr==maxyr,1,'last');
ERAi = find(yr==minyr,1,'first');
ERAf = find(yr==maxyr,1,'last');

tL = ERAf-ERAi+1;

% Do Only El Nino or only La Nina:
% $$$ n34(n34>-0.5) = 0;

% $$$ % Do Only IPO +ve or IPO -ve:
% $$$ IPODATA = load('tpi.timeseries.ersstv5.filt.dataraw');
% $$$ iponyrs = length(IPODATA(:,1));
% $$$ ipo = reshape(IPODATA(:,2:end)',[iponyrs*12 1]);
% $$$ ipo(abs(ipo)>2) = NaN;
% $$$ ipoyr = reshape(repmat(IPODATA(:,1),[1 12])',[iponyrs*12 1]);
% $$$ ipomn = repmat([1:12]',[iponyrs 1]);
% $$$ ipoi = find(ipoyr==minyr,1,'first');
% $$$ ipof = find(ipoyr==maxyr,1,'last');
% $$$ 
% $$$ inds = ipo(ipoi:ipof) > 0;
% $$$ 
% $$$ n34(n34i:n34f) = n34(n34i:n34f).*inds;

% Full year:
U10reg = reshape(reshape(U10(:,:,ERAi:ERAf),[xL*yL tL])*(n34(n34i: ...
                                                  n34f)-mean(n34(n34i:n34f)))/tL/std(n34(n34i:n34f)),[xL yL]);
V10reg = reshape(reshape(V10(:,:,ERAi:ERAf),[xL*yL tL])*(n34(n34i: ...
                                                  n34f)-mean(n34(n34i:n34f)))/tL/std(n34(n34i:n34f)),[xL yL]);
MSLreg = reshape(reshape(MSL(:,:,ERAi:ERAf),[xL*yL tL])*(n34(n34i: ...
                                                  n34f)-mean(n34(n34i:n34f)))/tL/std(n34(n34i:n34f)),[xL yL]);

%%%% Plotting Annual mean Nino 3.4:
figure;
% $$$ contourf(X,Y,U10reg,[-5:0.1:5],'linestyle','none');
contourf(X,Y,MSLreg,[-100 -5:0.25:5 100],'linestyle','none');
hold on;
nx = 6;ny = 6;
quiver(X(1:nx:xL,1:ny:yL),Y(1:nx:xL,1:ny:yL),U10reg(1:nx:xL,1:ny:yL),V10reg(1:nx:xL,1:ny:yL),1.5,'color','k');
plot(coastlon,coastlat,'-k')
% $$$ contour(X,Y,MSLreg,[25:25:500],'-k');
% $$$ contour(X,Y,MSLreg,[-500:25:-25],'--k');
caxis([-2.5 2.5]);
colormap(redblue)
cb = colorbar;
ylabel(cb,'mbar');
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
title([label ' Nino 3.4 (OISST 1982-2016) Regressed SLP and Winds']);
set(gca,'FontSize',25);

%%%% Month-by-Month:
U10regM = zeros(xL,yL,12);V10regM = U10regM; MSLregM = U10regM;
for i=1:12
    inds = find(mn(ERAi:ERAf) == i);
    n34inds = find(n34mn(n34i:n34f) == i);
    U10regM(:,:,i) = reshape(reshape(U10(:,:,inds+ERAi-1),[xL*yL length(inds)])*(n34(n34inds+n34i-1)-mean(n34(n34inds+n34i-1)))/length(inds)/std(n34(n34inds+n34i-1)),[xL yL]);
    V10regM(:,:,i) = reshape(reshape(V10(:,:,inds+ERAi-1),[xL*yL length(inds)])*(n34(n34inds+n34i-1)-mean(n34(n34inds+n34i-1)))/length(inds)/std(n34(n34inds+n34i-1)),[xL yL]);
    MSLregM(:,:,i) = reshape(reshape(MSL(:,:,inds+ERAi-1),[xL*yL length(inds)])*(n34(n34inds+n34i-1)-mean(n34(n34inds+n34i-1)))/length(inds)/std(n34(n34inds+n34i-1)),[xL yL]);
end


% Month-by-month:
figure;
for i=1:12
    subplot(4,3,i)
    contourf(X,Y,MSLregM(:,:,i),[-100 -5:0.25:5 100],'linestyle','none');
    hold on;
    nx = 10;ny = 10;
    quiver(X(1:nx:xL,1:ny:yL),Y(1:nx:xL,1:ny:yL),U10regM(1:nx:xL,1:ny:yL,i),V10regM(1:nx:xL,1:ny:yL,i),1.5,'color','k');
    plot(coastlon,coastlat,'-k')
% $$$ contour(X,Y,MSLregM(:,:,2*i),[25:25:500],'-k');
% $$$ contour(X,Y,MSLregM(:,:,2*i),[-500:25:-25],'--k');
    caxis([-2.5 2.5]);
    ylim([-90 -30]);
    colormap(redblue);
    cb = colorbar;
    ylabel(cb,'mbar');
% $$$     xlabel('Longitude ($^\circ$E)');
% $$$     ylabel('Latitude ($^\circ$N)');
    title(sprintf(['Nino 3.4 Reg' ...
                   'Month %02d'],i));
end

%%%% Particular Sequence:

figure;
inds = find((yr == 2015 & mn >= 7 ) | (yr == 2016 & mn <=6));
% $$$ inds = find((yr == 1997 & mn >= 7 ) | (yr == 1998 & mn <=6));
% $$$ inds = find((yr == 1982 & mn >= 7 ) | (yr == 1983 & mn <=6));
for i=1:length(inds)
    subplot(4,3,i)
    contourf(X,Y,MSL(:,:,inds(i)),[-100 -25:0.25:25 100],'linestyle','none');
    hold on;
    nx = 10;ny = 10;
    quiver(X(1:nx:xL,1:ny:yL),Y(1:nx:xL,1:ny:yL),U10(1:nx:xL,1:ny:yL,inds(i)),V10(1:nx:xL,1:ny:yL,inds(i)),1.5,'color','k');
    plot(coastlon,coastlat,'-k')
% $$$ contour(X,Y,MSLregM(:,:,2*i),[25:25:500],'-k');
% $$$ contour(X,Y,MSLregM(:,:,2*i),[-500:25:-25],'--k');
    caxis([-20 20]);
    ylim([-90 -30]);
    colormap(redblue);
    cb = colorbar;
    ylabel(cb,'mbar');
% $$$     xlabel('Longitude ($^\circ$E)');
% $$$     ylabel('Latitude ($^\circ$N)');
    title(sprintf(['%02d-%04d'],mn(inds(i)),yr(inds(i))));
end

%%% ASL -IPO - Nino 3.4 correlations:
figure;
% $$$ plot(yrd,ASL-mean(ASL),'-k');
plot(yrd,filter_field(ASL-mean(ASL),5*12+1,'-t'),'-k','linewidth',2);
hold on;
plot(yrd,filter_field(ASL-mean(ASL),10*12+1,'-t'),'--k','linewidth',2);
plot(yrd,filter_field(ASL-mean(ASL),15*12+1,'-t'),':k','linewidth',2);
hold on;
plot(n34yrd,filter_field(n34,5*12+1,'-t'),'-r','linewidth',2);
plot(ipoyrd,ipo,'-b','linewidth',2);
xlim([1980 2018]);
ylim([-1 1]);
legend('5-year low-pass ASL (mbar)','10-year low-pass ASL (mbar)','15-year low-pass ASL (mbar)','5-year low-pass Nino 3.4','IPO (TPI)');
xlabel('Year');
ylabel('SST ($^\circ$C) or mbar');

yr1 = 15;yr2 = 31; %1994-2010
yr1 = 13;yr2 = 32; %1992-2011
yr1 = 11;yr2 = 35; %1990-2014
yr1 = 19;yr2 = 36; %1998-2015
hold on;
plot([yrs(yr1) yrs(yr2)],[-0.2 -0.2-(yrs(yr1)-yrs(yr2))*trends(yr1,yr2)],'-y','linewidth',3);
text((yrs(yr1)+yrs(yr2))/2,-0.8,[num2str(yrs(yr1)) '-' num2str(yrs(yr2))],'color','y')

%%% ASL maximum trends:
yrs = unique(yr);
yrL = length(yrs);
[Y1,Y2] = ndgrid(yrs,yrs);

trends = NaN*zeros(size(Y1));
diffs = NaN*zeros(size(Y1));
for i=1:yrL
    for j=(i+5):yrL
        inds = find( (yr >= yrs(i)) & (yr <=yrs(j)));
        trends(i,j) = ((dnum(inds)-mean(dnum(inds))) \ ASL(inds)')*365;
        diffs(i,j) = ((dnum(inds)-mean(dnum(inds))) \ ASL(inds)')*365*(j-i+1);
    end
end

figure;
% $$$ pcolPlot(Y1,Y2,diffs*100);
pcolPlot(Y1,Y2,trends*100);
hold on;
plot(yrs,yrs+5,'--k');
plot(yrs,yrs+10,'--k');
plot(yrs,yrs+15,'--k');
plot(yrs,yrs+20,'--k');
plot(yrs,yrs+25,'--k');
plot(yrs,yrs+30,'--k');
caxis([-10 10]);
cb = colorbar;
ylabel(cb,'Pa/year');
xlabel('Initial year');
ylabel('Final year');
xlim([1980 2012]);
ylim([1985 2017]);

trends15 = trends;
for i=1:yrL
    if ((i+14)<=yrL)
        trends15(i,1:(i+14)) = NaN;
    end
end

% 20 year trends:
tr20 = zeros(yrL,1);
for i=1:(yrL-19)
    tr20(i) = trends(i,i+19);
end
figure;
plot(yrs,tr20*100,'-xk');
xlabel('Initial year');
ylabel('20-year ASL trend (Pa/yr)');

% 17 year trends:
tr17 = zeros(yrL,1);
for i=1:(yrL-19)
    tr17(i) = trends(i,i+16);
end
figure;
plot(yrs,tr17*100,'-xk');
xlabel('Initial year');
ylabel('17-year ASL trend (Pa/yr)');


%%% Long-term trend:
% $$$ mintr = zeros(1,1);
% $$$ for yr1=1980:1997
yr1 = 1991;
yr2 = 2010;
% $$$ yr2 = yr1+20-1;
inds = find( (yr >= yr1) & (yr <=yr2));

tLi = length(inds);

U10rs = reshape(double(U10(:,:,inds)),[xL*yL tLi]);
V10rs = reshape(double(V10(:,:,inds)),[xL*yL tLi]);
MSLrs = reshape(double(MSL(:,:,inds)),[xL*yL tLi]);
U10tr = reshape(((dnum(inds)-mean(dnum(inds))) \ (U10rs'))',[xL yL])*365;
V10tr = reshape(((dnum(inds)-mean(dnum(inds))) \ (V10rs'))',[xL yL])*365;
MSLtr = reshape(((dnum(inds)-mean(dnum(inds))) \ (MSLrs'))',[xL yL])*365;

% $$$ [tmp ln1] = min(abs(lon-200));
% $$$ [tmp ln2] = min(abs(lon-300));
% $$$ [tmp lt1] = min(abs(lat+80));
% $$$ [tmp lt2] = min(abs(lat+45));
% $$$ mintr(yr1-1980+1) = min(min(MSLtr(ln1:ln2,lt1:lt2)*100));
% $$$ [num2str(yr1) '-' num2str(yr2) ' ' num2str(min(min(MSLtr*100)))]
% $$$ end
% $$$ 
% $$$ hold on;
% $$$ plot(1980:1997,mintr,'-xr');
% $$$ xlabel('Initial year');
% $$$ ylabel('Minimum 20-year trend (Pa/yr) in ASL area');


figure;
contourf(X,Y,MSLtr*100,[-500 -100:2.5:100 500],'linestyle','none');
hold on;
nx = 6;ny = 6;
quiver(X(1:nx:xL,1:ny:yL),Y(1:nx:xL,1:ny:yL),U10tr(1:nx:xL,1:ny:yL),V10tr(1:nx:xL,1:ny:yL),1.5,'color','k');
quiver(X(1:nx:xL,1:ny:yL),Y(1:nx:xL,1:ny:yL),U10tr(1:nx:xL,1:ny:yL),V10tr(1:nx:xL,1:ny:yL),1.5,'color','k');
plot(coastlon,coastlat,'-k')
caxis([-30 30]);
colormap(redblue)
cb = colorbar;
ylabel(cb,'Pa yr^{-1}');
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
title([label ' SLP (Pa/year) and $10$m wind trends ' num2str(yr1) ...
       '-' num2str(yr2)]);
set(gca,'FontSize',25);

