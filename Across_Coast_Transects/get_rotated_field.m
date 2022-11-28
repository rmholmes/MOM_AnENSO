function [cc_field, lc_field] = get_rotated_field(fname,xvar,yvar,ti,lonr,latr,gridtype,vtype,angle)
% Gets a field on the new rotated grid.
dmn_bdy = 0;
[xL,yL] = size(lonr);

%Time average first:
if (length(ti) == 2)
    tL = ti(2)-ti(1)+1;
    ti = ti(1);
else
    tL = 1;
end

if (strfind(fname,'sose'))
    zL = 42; % max depth
    lonv = ncread(fname,'longitude')-360;
    latv = ncread(fname,'latitude'); %SOSE
else
    zL = 50;
    lonv = ncread(fname,'xt_ocean');
    latv = ncread(fname,'yt_ocean');
end
[lon,lat] = ndgrid(lonv,latv);

%%SOSE:

[tmp mnln] = min(abs(lonv-min(min(lonr))+dmn_bdy));
[tmp mxln] = min(abs(lonv-max(max(lonr))-dmn_bdy));
[tmp mnlt] = min(abs(latv-min(min(latr))+dmn_bdy));
[tmp mxlt] = min(abs(latv-max(max(latr))-dmn_bdy));
xiL = mxln-mnln+1;
etaL = mxlt-mnlt+1; %small domain lengths

lon = lon(mnln:mxln,mnlt:mxlt);
lat = lat(mnln:mxln,mnlt:mxlt);%shrink lon/lat.

if (vtype == '3D')
    xv = mean(ncread(fname,xvar,[mnln mnlt 1 ti],[xiL etaL zL tL]),4);
    if (~strcmp(yvar,''))
        yv = mean(ncread(fname,yvar,[mnln mnlt 1 ti],[xiL etaL zL tL]),4);
        lc_field = zeros(xL,yL,zL);
        for i=1:zL
            lc_field(:,:,i) = interp2(lon',lat',yv(:,:,i)',lonr,latr,'linear');
        end
    end

    cc_field = zeros(xL,yL,zL);
    for i=1:zL
        cc_field(:,:,i) = interp2(lon',lat',xv(:,:,i)',lonr,latr,'linear');
    end

elseif (vtype == '2D')
    xv = mean(ncread(fname,xvar,[mnln mnlt ti],[xiL etaL tL]),3);
    if (~strcmp(yvar,''))
        yv = mean(ncread(fname,yvar,[mnln mnlt ti],[xiL etaL tL]),3);
        lc_field = interp2(lon',lat',yv',lonr,latr,'linear');
    end
    cc_field = interp2(lon',lat',xv',lonr,latr,'linear');

elseif (vtype == 'h')
    xv = ncread(fname,xvar,[mnln mnlt],[xiL etaL]);
    cc_field = interp2(lon',lat',xv',lonr,latr,'linear');
end

if (~strcmp(yvar,''))
    tmp = cc_field;
    cc_field = sin(atan2(lc_field,tmp)-angle).*sqrt(tmp.^2+lc_field.^2);
    lc_field = cos(atan2(lc_field,tmp)-angle).*sqrt(tmp.^2+lc_field.^2);
else
    lc_field = 0;
end





