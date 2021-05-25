
% This script makes monthly average JRA55 netcdf files

clear all;
close all;

addpath(genpath('/short/e14/rmh561/software/matlab-utilities/'));
startup;

base = '/g/data/ua8/JRA55-do/v1-3/';
baseout = '/short/e14/rmh561/MOM_AnENSO/JRA55/';

U10in = [base sprintf('u_10.%04d.18Oct2017.nc',1980)];
V10in = [base sprintf('v_10.%04d.18Oct2017.nc',1980)];
slpin = [base sprintf('slp.%04d.18Aug2017.nc',1980)];

U10i = nc_inq(U10in,0);
U10i.Dimensions(3).Length = 1;
U10i.Variables(3).Dimensions(1).Length = 1;
U10i.Variables(3).Dimensions(1).Size = 1;
U10i.Variables(5).Dimensions(3).Length = 1;
U10i.Variables(5).Dimensions(3).Size(3) = 1;
U10i.Variables(5).Size(3) = 1;

V10i = nc_inq(V10in,0);
V10i.Dimensions(3).Length = 1;
V10i.Variables(3).Dimensions(1).Length = 1;
V10i.Variables(3).Dimensions(1).Size = 1;
V10i.Variables(5).Dimensions(3).Length = 1;
V10i.Variables(5).Dimensions(3).Size(3) = 1;
V10i.Variables(5).Size(3) = 1;

slpi = nc_inq(slpin,0);
slpi.Dimensions(3).Length = 1;
slpi.Variables(3).Dimensions(1).Length = 1;
slpi.Variables(3).Dimensions(1).Size = 1;
slpi.Variables(5).Dimensions(3).Length = 1;
slpi.Variables(5).Dimensions(3).Size(3) = 1;
slpi.Variables(5).Size(3) = 1;

lon = ncread(U10in,'longitude');
xL = length(lon);
lat = ncread(U10in,'latitude');
yL = length(lat);

for yr=1980:2016
    
    U10in = [base sprintf('u_10.%04d.18Oct2017.nc',yr)];
    V10in = [base sprintf('v_10.%04d.18Oct2017.nc',yr)];
    slpin = [base sprintf('slp.%04d.18Aug2017.nc',yr)];

    U10time = ncread(U10in,'time');
    U10timev = datevec(U10time+datenum([1900 1 1 0 0 0]));

    V10time = ncread(V10in,'time');
    V10timev = datevec(V10time+datenum([1900 1 1 0 0 0]));

    slptime = ncread(slpin,'time');
    slptimev = datevec(slptime+datenum([1900 1 1 0 0 0]));

    for mn=1:12
        sprintf('Doing Year %04d Month %02d',yr,mn)

        ind1 = find(U10timev(:,2) == mn,1,'first');
        ind2 = find(U10timev(:,2) == mn,1,'last');
        U10 = mean(ncread(U10in,'uas_10m',[1 1 ind1],[xL yL ind2-ind1+1]),3);
        U10i.FileName = [baseout sprintf('U10_%04d_%02d.nc',yr,mn)];
        ncid = nc_create([baseout sprintf('U10_%04d_%02d.nc',yr,mn)],'64bit_offset',U10i);
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'longitude'),lon);
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'latitude'),lat);
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'time'),[0],[1], ...
                      mean(U10time(ind1:ind2)));
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'uas_10m'),U10);
        netcdf.close(ncid);
        
        ind1 = find(V10timev(:,2) == mn,1,'first');
        ind2 = find(V10timev(:,2) == mn,1,'last');
        V10 = mean(ncread(V10in,'vas_10m',[1 1 ind1],[xL yL ind2-ind1+1]),3);
        V10i.FileName = [baseout sprintf('V10_%04d_%02d.nc',yr,mn)];
        ncid = nc_create([baseout sprintf('V10_%04d_%02d.nc',yr,mn)],'64bit_offset',V10i);
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'longitude'),lon);
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'latitude'),lat);
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'time'),[0],[1], ...
                      mean(V10time(ind1:ind2)));
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'vas_10m'),V10);
        netcdf.close(ncid);

        ind1 = find(slptimev(:,2) == mn,1,'first');
        ind2 = find(slptimev(:,2) == mn,1,'last');
        slp = mean(ncread(slpin,'psl',[1 1 ind1],[xL yL ind2-ind1+1]),3);
        slpi.FileName = [baseout sprintf('slp_%04d_%02d.nc',yr,mn)];
        ncid = nc_create([baseout sprintf('slp_%04d_%02d.nc',yr,mn)],'64bit_offset',slpi);
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'longitude'),lon);
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'latitude'),lat);
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'time'),[0],[1], ...
                      mean(slptime(ind1:ind2)));
        netcdf.putVar(ncid,netcdf.inqVarID(ncid,'psl'),slp);
        netcdf.close(ncid);

    end
end

    
