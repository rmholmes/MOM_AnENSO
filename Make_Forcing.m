
% This script makes forcing files for MOM-AnENSO runs.

clear all;
close all;

addpath(genpath('/g/data/e14/rmh561/software/matlab-utilities/'));
startup;

% Source data:
source_type = 1; % 0 = ERA, 1 = JRA
if (source_type)
    base = 'JRAdata/';
    label = 'JRA55-do';
else
    base = 'ERAdata/';
    label = 'ERA Interim';
end

% Output data:
out_type = 1; % 0 = ERA, 1 = JRA

%%% Get SST indices:

% $$$ % Nino 3.4:
% $$$ % OISST:
% $$$ DATA = load('index_data/NinoIndices.txt');
% $$$ n34 = DATA(:,10);
% $$$ n34yr = DATA(:,1);
% $$$ n34mn = DATA(:,2);
% $$$ n34yrd = n34yr + n34mn/12;

% $$$ % ERSST:
% $$$ DATA = load('index_data/nino34.data_ERSSTraw');
% HadISST:
DATA = load('index_data/nino34.long.data_HadISSTraw');
n34nyrs = length(DATA(:,1));
n34 = reshape(DATA(:,2:end)',[n34nyrs*12 1]);
n34(abs(n34)>50) = NaN;
n34yr = reshape(repmat(DATA(:,1),[1 12])',[n34nyrs*12 1]);
n34mn = repmat([1:12]',[n34nyrs 1]);
n34yrd = n34yr + n34mn/12;
n34cli = zeros(12,1);
for mi = 1:12
    n34cli(mi) = nanmean(n34(n34mn == mi & (n34yr >= 1981 & n34yr <=2010)));
    n34(n34mn == mi) = n34(n34mn == mi) - n34cli(mi);
end

%%% Get wind data:

Uname = [base 'U10_anom.nc'];
Vname = [base 'V10_anom.nc'];
if (source_type)
    U10 = ncread(Uname,'uas_10m');
    V10 = ncread(Vname,'vas_10m');
    lat = ncread(Uname,'latitude');
    lon = ncread(Uname,'longitude');
    time = ncread(Uname,'time');
    dnum = datenum([1900 1 1 0 0 0])+time;
else
    U10 = ncread(Uname,'10U_GDS0_SFC');
    V10 = ncread(Vname,'10V_GDS0_SFC');
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

%%% Calculate N34 regression on source data grid:
minyr = max([min(yr) min(n34yr)]);
maxyr = min([max(yr) max(n34yr)]);
minyr = 1982;

n34i = find(n34yr==minyr,1,'first');
n34f = find(n34yr==maxyr,1,'last');
ERAi = find(yr==minyr,1,'first');
ERAf = find(yr==maxyr,1,'last');

tL = ERAf-ERAi+1;

% Full year:
U10reg = reshape(reshape(U10(:,:,ERAi:ERAf),[xL*yL tL])*(n34(n34i: ...
                                                  n34f)-mean(n34(n34i:n34f)))/tL/std(n34(n34i:n34f)),[xL yL]);
V10reg = reshape(reshape(V10(:,:,ERAi:ERAf),[xL*yL tL])*(n34(n34i: ...
                                                  n34f)-mean(n34(n34i:n34f)))/tL/std(n34(n34i:n34f)),[xL yL]);

% $$$ %%% Specified period trends:
% $$$ yr1 = 1991;
% $$$ yr2 = 2010;
% $$$ inds = find( (yr >= yr1) & (yr <=yr2));
% $$$ tLi = length(inds);
% $$$ U10rs = reshape(double(U10(:,:,inds)),[xL*yL tLi]);
% $$$ V10rs = reshape(double(V10(:,:,inds)),[xL*yL tLi]);
% $$$ U10tr = reshape(((dnum(inds)-mean(dnum(inds))) \ (U10rs'))',[xL yL]);
% $$$ V10tr = reshape(((dnum(inds)-mean(dnum(inds))) \ (V10rs'))',[xL yL]);
% $$$ 
% $$$ CNYFtime = ncread('CNYFu_10.nc','TIME');
% $$$ tL = length(CNYFtime);
% $$$ 
% $$$ CNYFlon = ncread('CNYFu_10.nc','LON');
% $$$ xL = length(CNYFlon);
% $$$ CNYFlat = ncread('CNYFu_10.nc','LAT');
% $$$ yL = length(CNYFlat);
% $$$ [CX,CY] = ndgrid(CNYFlon,CNYFlat);
% $$$ 
% $$$ U10trCNYF = interp2(X',Y',U10tr',CX,CY,'linear');
% $$$ V10trCNYF = interp2(X',Y',V10tr',CX,CY,'linear');
% $$$ 
% $$$ CNYFu10 = ncread('CNYFu_10.nc','U_10_MOD');
% $$$ CNYFv10 = ncread('CNYFv_10.nc','V_10_MOD');
% $$$ 
% $$$ outfold = 'ASL_20yrtrend_global/';
% $$$ 
% $$$ for yi = 1:20
% $$$     sprintf('Doing %02d',yi)
% $$$     outnameU = [outfold 'CNYFu_10_ASL20yrglobal_yr' ...
% $$$                         num2str(yi) '.nc'];
% $$$     outnameV = [outfold 'CNYFv_10_ASL20yrglobal_yr' ...
% $$$                         num2str(yi) '.nc'];
% $$$     copyfile('CNYFu_10.nc',outnameU);
% $$$     copyfile('CNYFv_10.nc',outnameV);
% $$$     
% $$$     U10wt = CNYFu10 + repmat(U10trCNYF,[1 1 tL]).*repmat(permute(CNYFtime+365*(yi-1),[2 ...
% $$$                         3 1]),[xL yL 1]);
% $$$     V10wt = CNYFv10 + repmat(V10trCNYF,[1 1 tL]).*repmat(permute(CNYFtime+365*(yi-1),[2 ...
% $$$                         3 1]),[xL yL 1]);
% $$$ 
% $$$     ncid = netcdf.open(outnameU,'NC_WRITE');
% $$$     netcdf.putVar(ncid,netcdf.inqVarID(ncid,'U_10_MOD'),U10wt);
% $$$     netcdf.close(ncid);
% $$$ 
% $$$     ncid = netcdf.open(outnameV,'NC_WRITE');
% $$$     netcdf.putVar(ncid,netcdf.inqVarID(ncid,'V_10_MOD'),V10wt);
% $$$     netcdf.close(ncid);
% $$$ end
% $$$ 
% $$$ % Diff for check:
% $$$ for yi = 1:20
% $$$     outnameU = [outfold 'CNYFu_10_ASL20yrglobal_yr' ...
% $$$                         num2str(yi) '.nc'];
% $$$     outnameV = [outfold 'CNYFv_10_ASL20yrglobal_yr' ...
% $$$                         num2str(yi) '.nc'];
% $$$     system(['ncdiff ' outnameU ' CNYFu_10.nc ' outnameU(1:end-3) ...
% $$$             '_diff.nc']);
% $$$     system(['ncdiff ' outnameV ' CNYFv_10.nc ' outnameV(1:end-3) ...
% $$$             '_diff.nc']);
% $$$ end
% $$$ 


% $$$ %%%%%%%%%%%%% Realistic N34 time series:
% $$$ yr1 = 1968;
% $$$ yr2 = 2017;
% $$$ CNYFtime = ncread('CNYFu_10.nc','TIME')+datenum([1900 1 1 0 0 0]);
% $$$ tL = length(CNYFtime);
% $$$ tCNYF = linspace(0,12,tL);
% $$$ tMON = -0.5:1:12.5;
% $$$ 
% $$$ CNYFlon = ncread('CNYFu_10.nc','LON');
% $$$ xL = length(CNYFlon);
% $$$ CNYFlat = ncread('CNYFu_10.nc','LAT');
% $$$ yL = length(CNYFlat);
% $$$ [CX,CY] = ndgrid(CNYFlon,CNYFlat);
% $$$ 
% $$$ U10regCNYF = interp2(X',Y',U10reg',CX,CY,'linear');
% $$$ V10regCNYF = interp2(X',Y',V10reg',CX,CY,'linear');
% $$$ 
% $$$ CNYFu10 = ncread('CNYFu_10.nc','U_10_MOD');
% $$$ CNYFv10 = ncread('CNYFv_10.nc','V_10_MOD');
% $$$ 
% $$$ % Copied from /short/v45/pas561/mom/input/core2iaf/u_10.1948-2007.06JUN2011.nc
% $$$ % Copied from /short/v45/pas561/mom/input/core2iaf/v_10.1948-2007.06JUN2011.nc
% $$$ outnameU = 'IIAF/u_10_iaf_n34ideal.nc';
% $$$ outnameV = 'IIAF/v_10_iaf_n34ideal.nc';
% $$$ 
% $$$ for yi = yr1:yr2
% $$$     sprintf('Doing %02d',yi)
% $$$     tini = find(n34yr == yi,1,'first');
% $$$     tfin = find(n34yr == yi,1,'last');
% $$$     N34 = interp1(tMON,[n34(tini-1:tfin+1)],tCNYF);
% $$$     
% $$$     U10 = CNYFu10 + repmat(U10regCNYF,[1 1 tL]).*repmat(permute(N34,[1 ...
% $$$                         3 2]),[xL yL 1]);
% $$$     V10 = CNYFv10 + repmat(V10regCNYF,[1 1 tL]).*repmat(permute(N34,[1 ...
% $$$                         3 2]),[xL yL 1]);
% $$$ 
% $$$     ncid = netcdf.open(outnameU,'NC_WRITE');
% $$$     netcdf.putVar(ncid,netcdf.inqVarID(ncid,'U_10_MOD'),[0 0 ...
% $$$                         (yi-yr1)*1460],[xL yL 1460],U10);
% $$$     netcdf.close(ncid);
% $$$ 
% $$$     ncid = netcdf.open(outnameV,'NC_WRITE');
% $$$     netcdf.putVar(ncid,netcdf.inqVarID(ncid,'V_10_MOD'),[0 0 ...
% $$$                         (yi-yr1)*1460],[xL yL 1460],V10);
% $$$     netcdf.close(ncid);
% $$$ end
% $$$ %Zero out other years:
% $$$ tL = length(ncread(outnameU,'TIME'));
% $$$ ncid = netcdf.open(outnameU,'NC_WRITE');
% $$$ netcdf.putVar(ncid,netcdf.inqVarID(ncid,'U_10_MOD'),[0 0 ...
% $$$                     (yr2+1-yr1)*1460],[xL yL tL-(yr2+1-yr1)*1460],zeros(xL,yL,tL-(yr2+1-yr1)*1460));
% $$$ netcdf.close(ncid);
% $$$ 
% $$$ ncid = netcdf.open(outnameV,'NC_WRITE');
% $$$ netcdf.putVar(ncid,netcdf.inqVarID(ncid,'V_10_MOD'),[0 0 ...
% $$$                     (yr2+1-yr1)*1460],[xL yL tL-(yr2+1-yr1)*1460],zeros(xL,yL,tL-(yr2+1-yr1)*1460));
% $$$ netcdf.close(ncid);
% $$$ 
% $$$ % Check one point:
% $$$ [tmp lnpt] = min(abs(CNYFlon-170));
% $$$ [tmp ltpt] = min(abs(CNYFlat));
% $$$ 
% $$$ winds = squeeze(ncread(outnameU,'U_10_MOD',[lnpt ltpt 1],[1 1 tL]));
% $$$ windsbase = squeeze(CNYFu10(lnpt,ltpt,:));
% $$$ anoms = winds - repmat(windsbase,[60 1]);
% $$$ 
% $$$ figure;
% $$$ plot((1:length(anoms))/(4*365),anoms,'-k');
% $$$ hold on;
% $$$ plot(n34yrd(n34yr>=yr1 & n34yr <= yr2)-yr1,n34(n34yr>=yr1 & n34yr <= yr2)*1.7,'-r');


%%%%%%%%%%%%% Idealized experiments:
% Plot time series and fit:
% $$$ 
% Extreme El Nino's:
inds = find((n34yr == 2015) | (n34yr == 2016) | (n34yr == 2017));
n1516 = n34(inds);
inds = find((n34yr == 1997) | (n34yr == 1998) | (n34yr == 1999));
n9798 = n34(inds);
inds = find((n34yr == 1982) | (n34yr == 1983) | (n34yr == 1984));
n8283 = n34(inds);

avg = (n8283+n1516+n9798)/3;
outname = '_ELNINO';

% $$$ % Extreme La Nina's:
% $$$ inds = find((n34yr == 1988) | (n34yr == 1989) | (n34yr == 1990));
% $$$ n8889 = n34(inds);
% $$$ inds = find((n34yr == 2007) | (n34yr == 2008) | (n34yr == 2009));
% $$$ n0708 = n34(inds);
% $$$ inds = find((n34yr == 2010) | (n34yr == 2011) | (n34yr == 2012));
% $$$ n1011 = n34(inds);

% $$$ avg = (n8889+n0708+n1011)/3;
% $$$ outname = '_LANINA';

avgs = avg;
avgs(30:end) = 0;
if (strcmp(outname,'_LANINA'))
    avgs(1:10) = min(avgs(1:10),0);
end
hfilt = 1;
avgsext = [zeros(hfilt,1); avgs; zeros(hfilt,1)];
avgsext = filter_field(avgsext,2*hfilt+1,'-t');
avgs = avgsext(hfilt+1:end-hfilt);
avgs(1) = 0;


% $$$ clf;
% $$$ % $$$ figure;
% $$$ set(gcf,'Position',[2100         235        1087         545]);
% $$$ % $$$ plot(1:36,n1516,'-m','linewidth',2);
% $$$ % $$$ hold on;
% $$$ % $$$ plot(1:36,n9798,'-b','linewidth',2);
% $$$ % $$$ plot(1:36,n8283,'-r','linewidth',2);
% $$$ plot(1:36,n8889,'-m','linewidth',2);
% $$$ hold on;
% $$$ plot(1:36,n0708,'-b','linewidth',2);
% $$$ plot(1:36,n1011,'-r','linewidth',2);
% $$$ plot(1:36,avg,'-c','linewidth',2);
% $$$ plot(1:36,avgs,'-k','linewidth',4);
% $$$ % $$$ legend('2015-2017','1997-1999','1982-1984','Average','3-month Smoothed Average');
% $$$ legend('1988-1990','2007-2009','2010-2012','Average','3-month Smoothed Average');
% $$$ xlabel('Month');
% $$$ % $$$ set(gca,'xtick',[1:12]);
% $$$ % $$$ set(gca,'xticklabel',[7:12 1:6]);
% $$$ ylabel('Nino 3.4 Anomaly $(^\circ$C)');

% $$$ % Read from MOM:
% $$$ base = '/g/data/e14/rmh561/MOM_AnENSO/';
% $$$ lon = ncread([base 'temp.cat.diff_yr1.nc'],'xt_ocean');
% $$$ lat = ncread([base 'temp.cat.diff_yr1.nc'],'yt_ocean');
% $$$ [tmp ind1] = min(abs(lon+170));
% $$$ [tmp ind2] = min(abs(lon+120));
% $$$ [tmp ind3] = min(abs(lat+5));
% $$$ [tmp ind4] = min(abs(lat-5));
% $$$ 
% $$$ time = [ncread([base 'temp.cat.diff_yr1.nc'],'time'); ...
% $$$         ncread([base 'temp.cat.diff_yr2.nc'],'time'); ...
% $$$         ncread([base 'temp.cat.diff_yr3.nc'],'time')];
% $$$ tL = length(time);
% $$$ temp = [squeeze(nanmean(nanmean(ncread([base 'temp.cat.diff_yr1.nc'],'temp',[ind1 ind3 1 1],[ind2-ind1+1 ind4-ind3+1 1 tL/3]),1),2)); ...
% $$$         squeeze(nanmean(nanmean(ncread([base 'temp.cat.diff_yr2.nc'],'temp',[ind1 ind3 1 1],[ind2-ind1+1 ind4-ind3+1 1 tL/3]),1),2)); ...
% $$$         squeeze(nanmean(nanmean(ncread([base 'temp.cat.diff_yr3.nc'],'temp',[ind1 ind3 1 1],[ind2-ind1+1 ind4-ind3+1 1 tL/3]),1),2))];
% $$$ 
% $$$ time_mon = (time-time(1)+15)/365*12;
% $$$ hold on;
% $$$ plot(time_mon,temp,'--k','linewidth',2)
% $$$ legend('2015-2017','1997-1999','1982-1984','Average','3-month Smoothed Average','MOM025');


if (out_type == 0)
    % Make new CORE-NYF u10/v10:

% $$$ baseCNYF = '/short/e14/rmh561/mom/input/gfdl_nyf_1080_clean/';
% $$$ copyfile([baseCNYF 'u_10.nc'],'CNYFu_10.nc');
% $$$ copyfile([baseCNYF 'v_10.nc'],'CNYFv_10.nc');

    CNYFtime = ncread('CNYFu_10.nc','TIME')+datenum([1900 1 1 0 0 0]);
    tL = length(CNYFtime);
    tCNYF = linspace(0,12,tL);
    tMON = -0.5:1:12.5;

    CNYFlon = ncread('CNYFu_10.nc','LON');
    xL = length(CNYFlon);
    CNYFlat = ncread('CNYFu_10.nc','LAT');
    yL = length(CNYFlat);
    [CX,CY] = ndgrid(CNYFlon,CNYFlat);

    U10regCNYF = interp2(X',Y',U10reg',CX,CY,'linear');
    V10regCNYF = interp2(X',Y',V10reg',CX,CY,'linear');

    CNYFu10 = ncread('CNYFu_10.nc','U_10_MOD');
    CNYFv10 = ncread('CNYFv_10.nc','V_10_MOD');

% Idealized three-year runs:
for yi = 1:3
    sprintf('Doing %02d',yi)
    outnameU = sprintf(['u_10' outname '_yr%01d.nc'],yi);
    outnameV = sprintf(['v_10' outname '_yr%01d.nc'],yi);
    copyfile('CNYFu_10.nc',outnameU);
    copyfile('CNYFv_10.nc',outnameV);
    
    if yi == 1
        N34 = interp1(tMON,[0; avgs(1:13)],tCNYF);
    elseif yi == 2
        N34 = interp1(tMON,avgs(12:25),tCNYF);
    elseif yi == 3
        N34 = interp1(tMON,[avgs(24:36); 0],tCNYF);
    end         
    
    U10 = CNYFu10 + repmat(U10regCNYF,[1 1 tL]).*repmat(permute(N34,[1 ...
                        3 2]),[xL yL 1]);
    V10 = CNYFv10 + repmat(V10regCNYF,[1 1 tL]).*repmat(permute(N34,[1 ...
                        3 2]),[xL yL 1]);

    ncid = netcdf.open(outnameU,'NC_WRITE');
    netcdf.putVar(ncid,netcdf.inqVarID(ncid,'U_10_MOD'),U10);
    netcdf.close(ncid);

    ncid = netcdf.open(outnameV,'NC_WRITE');
    netcdf.putVar(ncid,netcdf.inqVarID(ncid,'V_10_MOD'),V10);
    netcdf.close(ncid);
end

% Diff for check:
for yi = 1:3
    outnameU = sprintf(['u_10' outname '_yr%01d.nc'],yi);
    outnameV = sprintf(['v_10' outname '_yr%01d.nc'],yi);
    system(['ncdiff ' outnameU ' CNYFu_10.nc ' outnameU(1:end-3) ...
            '_diff.nc']);
    system(['ncdiff ' outnameV ' CNYFv_10.nc ' outnameV(1:end-3) ...
            '_diff.nc']);
end


% $$$ % Full nino 3.4 time series:
% $$$ for yi = 1982:2016
% $$$     sprintf('Doing %02d',yi)
% $$$     outnameU = sprintf(['IIAF/u_10' outname '_yr%04d.nc'],yi);
% $$$     outnameV = sprintf(['IIAF/v_10' outname '_yr%04d.nc'],yi);
% $$$     copyfile('CNYFu_10.nc',outnameU);
% $$$     copyfile('CNYFv_10.nc',outnameV);
% $$$     
% $$$     if yi == 1982
% $$$         tini = find(n34yr == yi,1,'first');
% $$$         tfin = find(n34yr == yi,1,'last');
% $$$         N34 = interp1(tMON,[0; n34(tini:tfin+1)],tCNYF);
% $$$     else
% $$$         tini = find(n34yr == yi,1,'first');
% $$$         tfin = find(n34yr == yi,1,'last');
% $$$         N34 = interp1(tMON,[n34(tini-1:tfin+1)],tCNYF);
% $$$     end        
% $$$     
% $$$     U10 = CNYFu10 + repmat(U10regCNYF,[1 1 tL]).*repmat(permute(N34,[1 ...
% $$$                         3 2]),[xL yL 1]);
% $$$     V10 = CNYFv10 + repmat(V10regCNYF,[1 1 tL]).*repmat(permute(N34,[1 ...
% $$$                         3 2]),[xL yL 1]);
% $$$ 
% $$$     ncid = netcdf.open(outnameU,'NC_WRITE');
% $$$     netcdf.putVar(ncid,netcdf.inqVarID(ncid,'U_10_MOD'),U10);
% $$$     netcdf.close(ncid);
% $$$ 
% $$$     ncid = netcdf.open(outnameV,'NC_WRITE');
% $$$     netcdf.putVar(ncid,netcdf.inqVarID(ncid,'V_10_MOD'),V10);
% $$$     netcdf.close(ncid);
% $$$ end
% $$$ 
% $$$ % Diff for check:
% $$$ for yi = 1982:2016
% $$$     outnameU = sprintf(['IIAF/u_10' outname '_yr%04d.nc'],yi);
% $$$     outnameV = sprintf(['IIAF/v_10' outname '_yr%04d.nc'],yi);
% $$$     system(['ncdiff ' outnameU ' CNYFu_10.nc ' outnameU(1:end-3) ...
% $$$             '_diff.nc']);
% $$$     system(['ncdiff ' outnameV ' CNYFv_10.nc ' outnameV(1:end-3) ...
% $$$             '_diff.nc']);
% $$$ end


else
    % Make new JRA55 v1.3 u10/v10:

    baseJRA = '/g/data/ik11/inputs/JRA-55/RYF/v1-3/';
    u10fn = 'RYF.u_10.1990_1991.nc';
    v10fn = 'RYF.v_10.1990_1991.nc';
    copyfile([baseJRA u10fn],u10fn);
    copyfile([baseJRA v10fn],v10fn);

    JRAtime = ncread(u10fn,'time')+datenum([1900 1 1 0 0 0]);
    tL = length(JRAtime);
    tJRA = linspace(0,12,tL);
    tMON = -0.5:1:12.5;

% $$$     % No need for interpolation if using JRA55 source data:
% $$$     JRAlon = ncread(u10fn,'lon');
% $$$     xL = length(JRAlon);
% $$$     JRAlat = ncread(u10fn,'lat');
% $$$     yL = length(JRAlat);
% $$$     [CX,CY] = ndgrid(JRAlon,JRAlat);
% $$$ 
% $$$     U10regJRA = interp2(X',Y',U10reg',CX,CY,'linear');
% $$$     V10regJRA = interp2(X',Y',V10reg',CX,CY,'linear');
% $$$ 
    U10regJRA = U10reg;
    V10regJRA = V10reg;

    JRAu10 = ncread(u10fn,'uas_10m');
    JRAv10 = ncread(v10fn,'vas_10m');

% Idealized three-year runs:
for yi = 1:3
    sprintf('Doing %02d',yi)
    outnameU = sprintf(['RYF.u_10.1990_1991' outname '_yr%01d.nc'],yi);
    outnameV = sprintf(['RYF.v_10.1990_1991' outname '_yr%01d.nc'],yi);
    copyfile(u10fn,outnameU);
    copyfile(v10fn,outnameV);
    
    if yi == 1
        N34 = interp1(tMON,[0; avgs(1:13)],tJRA);
    elseif yi == 2
        N34 = interp1(tMON,avgs(12:25),tJRA);
    elseif yi == 3
        N34 = interp1(tMON,[avgs(24:36); 0],tJRA);
    end         
    
    U10 = JRAu10 + repmat(U10regJRA,[1 1 tL]).*repmat(permute(N34,[1 ...
                        3 2]),[xL yL 1]);
    V10 = JRAv10 + repmat(V10regJRA,[1 1 tL]).*repmat(permute(N34,[1 ...
                        3 2]),[xL yL 1]);

    ncid = netcdf.open(outnameU,'NC_WRITE');
    netcdf.putVar(ncid,netcdf.inqVarID(ncid,'uas_10m'),U10);
    netcdf.close(ncid);

    ncid = netcdf.open(outnameV,'NC_WRITE');
    netcdf.putVar(ncid,netcdf.inqVarID(ncid,'vas_10m'),V10);
    netcdf.close(ncid);
end

% Diff for check:
for yi = 1:3
    outnameU = sprintf(['RYF.u_10.1990_1991' outname '_yr%01d.nc'],yi);
    outnameV = sprintf(['RYF.v_10.1990_1991' outname '_yr%01d.nc'],yi);
    system(['ncdiff ' outnameU ' ' u10fn ' ' outnameU(1:end-3) ...
            '_diff.nc']);
    system(['ncdiff ' outnameV ' ' v10fn ' ' outnameV(1:end-3) ...
            '_diff.nc']);
end
    end
    
