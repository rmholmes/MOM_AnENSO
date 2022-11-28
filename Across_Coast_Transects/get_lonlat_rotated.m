function [lon_rot,lat_rot,Corners,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cn)
% Gets the coordinates of the rotated lon and lat coordinates at a
% center point cn = lon, lat, angle to north.
%
% W = 

    Re = 6378000.0;
    lat_to_km = 2*pi*Re/360.0/1e3;

    dw = (W+Wm)/Nw;
    dl = L/Nl;
    theta = cn(3)+pi/2; %angle of front to eastward
    lon_rot = zeros(Nw+1,Nl+1); 
    lat_rot = lon_rot; %initialize lon and lat of rotated grid.
    for i=1:(Nw+1)
        for j=1:(Nl+1)
            xp = (i-1)*dw-Wm; %x-spacing from center point
            yp = (j-(Nl/2+1))*dl; %y-spacing from center point
            phi = atan2(yp,xp); %angle from center point
            leng = sqrt(xp^2+yp^2); %distance from center point

            lat_rot(i,j) = cn(2)+sin(phi+theta)*leng/lat_to_km; %lat of point
            lon_rot(i,j) = cn(1)+cos(phi+theta)*leng/lat_to_km/ ...
                cos(pi/180*(lat_rot(i,j)+cn(2))/2); %lon of point
        end
    end
    
    %Domain corners:
    Corners = [lon_rot(1,1) lat_rot(1,1); lon_rot(1,end) lat_rot(1,end);...
               lon_rot(end,end) lat_rot(end,end); lon_rot(end,1) lat_rot(end,1)];

    %cross-coast dimension:
    cc = (0:dw:(Nw*dw))-Wm;
    cc = (0:Nw)*dw-Wm;
    
    %along-coast dimension:
    lc = ((1:(Nl+1))-(Nl/2+1))*dl;
   
end

