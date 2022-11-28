function [c,dx,cseg,dc] = get_isobath_sections(x_rho,y_rho,mask,h,isobath,model,plotting,sp,isp)
% Gets a smooth isobath from either MOM5 or ROMS output for
% Antarctica. Then averages this isobath into sections of length
% sp.

if (~plotting)
    f = figure('visible','off');
end
Re = 6378000.0;
lat_to_km = 2*pi*Re/360.0/1e3;
contour(x_rho,y_rho,mask,[0.5 0.5],'-k','LineWidth',1);
hold on;
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
h(isnan(h)) = 0;
[craw,tmp] = contour(x_rho,y_rho,h,isobath*[1 1],'-b','LineWidth',1);

if (strcmp(model,'MOM'))
    if (isobath == 1000)
% MOM 1000m isobath exclusions:
craw(:,abs(craw(1,:))>361.0) = [];
craw(:,abs(craw(2,:))>91.0) = [];
c = craw;
filtpts = 51;
c = c(:,c(2,:)<-59.9);
c = c(:,c(2,:)>-76.6);
inds = c(2,:)>-68.5 & c(2,:)<-65.0 & c(1,:)>-200.0 & c(1,:) < -175.0;c(:,inds) = [];
inds = c(2,:)>-71.85 & c(2,:)<-50.0 & c(1,:)>-185.0 & c(1,:) < -178.4;c(:,inds) = [];
inds = c(2,:)>-64.95 & c(2,:)<-50.0 & c(1,:)>-247.0 & c(1,:) < -246.0;c(:,inds) = [];
inds = c(2,:)>-72.1 & c(2,:)<-50.0 & c(1,:)>-119.0 & c(1,:) < -118.0;c(:,inds) = [];
inds = c(2,:)>-72.4 & c(2,:)<-71.30 & c(1,:)>-104.5 & c(1,:) < -100.0;c(:,inds) = [];
inds = c(2,:)>-69.5 & c(2,:)<-50.0 & c(1,:)>-95 & c(1,:) < -85;c(:,inds) = [];
inds = c(2,:)>-63.0 & c(2,:)<-62.0 & c(1,:)>-60.0 & c(1,:) < -54.0;c(:,inds) = [];
inds = c(2,:)>-62.1 & c(2,:)<-61.45 & c(1,:)>-57.0 & c(1,:) < -54;c(:,inds) = [];
inds = c(2,:)>-62.5 & c(2,:)<-50.0 & c(1,:)>-51.1 & c(1,:) < -20.0;c(:,inds) = [];
inds = c(2,:)>-60.9 & c(2,:)<-60.1 & c(1,:)>-52.0 & c(1,:) < -49.5;c(:,inds) = [];
inds = c(2,:)>-66.5 & c(2,:)<-50.0 & c(1,:)>30.0 & c(1,:) < 40.0;c(:,inds) = [];
inds = c(2,:)>-71.153 & c(2,:)<-71.148 & c(1,:)>-106.135 & c(1,:) < -106.11;c(:,inds) = [];
npts = length(c(1,:));
    elseif (isobath == 2500)
        % MOM 2500m isobath exclusions:
craw(:,abs(craw(1,:))>361.0) = [];
craw(:,abs(craw(2,:))>91.0) = [];
c = craw;
filtpts = 51;
c = c(:,c(2,:)<-60.0);
c = c(:,c(2,:)>-76.6);
inds = c(2,:)>-64 & c(2,:)<-0.0 & c(1,:)>-5000.0 & c(1,:) < -272.0;c(:,inds) = [];
inds = c(2,:)>-66.1 & c(2,:)<-0.0 & c(1,:)>-190.0 & c(1,:) < -80.0;c(:,inds) = [];
inds = c(2,:)>-65 & c(2,:)<-0.0 & c(1,:)>-210.0 & c(1,:) < -80.0;c(:,inds) = [];
inds = c(2,:)>-70.1 & c(2,:)<-0.0 & c(1,:)>-185.0 & c(1,:) < -160.0;c(:,inds) = [];
inds = c(2,:)>-70.68 & c(2,:)<-0.0 & c(1,:)>-135.0 & c(1,:) < -110.0;c(:,inds) = [];
inds = c(2,:)>-71.58 & c(2,:)<-0.0 & c(1,:)>-124.5 & c(1,:) < -121.0;c(:,inds) = [];
inds = c(2,:)>-69.2 & c(2,:)<-68.2 & c(1,:)>-92.5 & c(1,:) < -89.0;c(:,inds) = [];
inds = c(2,:)>-64 & c(2,:)<-63.85 & c(1,:)>-253.4 & c(1,:) < -252.8;c(:,inds) = [];
inds = c(2,:)>-73.8 & c(2,:)<-73.5 & c(1,:)>-173.0 & c(1,:) < -170.0;c(:,inds) = [];
inds = c(2,:)>-71.2 & c(2,:)<-70.0 & c(1,:)>-188.0 & c(1,:) < -186.5;c(:,inds) = [];
inds = c(2,:)>-68.5 & c(2,:)<-66.0 & c(1,:)>-189.9 & c(1,:) < -185;c(:,inds) = [];
inds = c(2,:)>-69.1 & c(2,:)<-68.6 & c(1,:)>-191 & c(1,:) < -189.5;c(:,inds) = [];
inds = c(2,:)>-69.2 & c(2,:)<-67.8 & c(1,:)>-195.5 & c(1,:) < -192;c(:,inds) = [];
inds = c(2,:)>-66.8 & c(2,:)<-0.0 & c(1,:)>-200 & c(1,:) < -199;c(:,inds) = [];
inds = c(2,:)>-65.8 & c(2,:)<-0.0 & c(1,:)>-200 & c(1,:) < -198;c(:,inds) = [];
inds = c(2,:)>-65.6 & c(2,:)<-0.0 & c(1,:)>-192 & c(1,:) < -188;c(:,inds) = [];
inds = c(2,:)>-67 & c(2,:)<-66.4 & c(1,:)>-196.5 & c(1,:) < -195.3;c(:,inds) = [];
inds = c(2,:)>-65.4 & c(2,:)<-65.1 & c(1,:)>-70.5 & c(1,:) < -70;c(:,inds) = [];
inds = c(2,:)>-60.7 & c(2,:)<-59.8 & c(1,:)>-68 & c(1,:) < -56.5;c(:,inds) = [];
inds = c(2,:)>-60.8 & c(2,:)<-60.225 & c(1,:)>-52 & c(1,:) < -49.5;c(:,inds) = [];
inds = c(2,:)>-64.53 & c(2,:)<-0.0 & c(1,:)>-40.81 & c(1,:) < -0.0;c(:,inds) = [];
inds = c(2,:)>-60.7 & c(2,:)<-0.0 & c(1,:)>-42 & c(1,:) < -0.0;c(:,inds) = [];
inds = c(2,:)>-63.6 & c(2,:)<-62.655 & c(1,:)>-49.3 & c(1,:) < -0.0;c(:,inds) = [];
inds = c(2,:)>-68.2 & c(2,:)<-67.5 & c(1,:)>-53.5 & c(1,:) < -52.5;c(:,inds) = [];
inds = c(2,:)>-69.3 & c(2,:)<-68.9 & c(1,:)>3.5 & c(1,:) < 4.25;c(:,inds) = [];
inds = c(2,:)>-70.815 & c(2,:)<-70.6 & c(1,:)>-121 & c(1,:) < -119;c(:,inds) = [];

inds = c(2,:)>-1000 & c(2,:)<0 & c(1,:)>0 & c(1,:) < 1000;c(:,inds) = [];
npts = length(c(1,:));
    elseif (isobath == 100)
        % 50m (coast isobath):
        craw(:,abs(craw(1,:))>361.0) = [];
        craw(:,abs(craw(2,:))>91.0) = [];
        craw(:,craw(1,:)>-78) = [];
        craw(:,craw(1,:)<-170) = [];
        craw(:,craw(2,:)>-71) = [];
        c = craw;
        filtpts = 31;
        inds = c(2,:)>-76 & c(2,:)<-75.75 & c(1,:)>-150.2 & c(1,:) < -150.0;c(:,inds) = [];
        npts = length(c(1,:));
    end
else
% ROMS 1000m isobath exclusions:
c = craw;
filtpts = 13;
c = c(:,c(2,:)<-61.0);
inds = c(2,:)<-77.8 & c(1,:)>300.0;c(:,inds) = [];
inds = c(2,:)>-63.8 & c(1,:)>311.7;c(:,inds) = [];
inds = c(2,:)>-72 & c(1,:)>240.0 & c(1,:) < 244.0;c(:,inds) = [];
inds = c(2,:)>-64.8 & c(2,:)<-64.2 & c(1,:)>298.0 & c(1,:) < 298.8;c(:,inds) = [];
inds = c(2,:)>-68.8 & c(2,:)<-67.4 & c(1,:)>165.0 & c(1,:) < 175.0;c(:,inds) = [];
inds = c(2,:)>-72.54 & c(2,:)<-72.46 & c(1,:)>251.4 & c(1,:) < 251.8;c(:,inds) = [];
npts = length(c(1,:));
end

%Smooth isobath:
cf = zeros(2,npts+filtpts*2);
cf(:,(filtpts+1):(filtpts+1+npts-1)) = c;
cf(:,1:filtpts) = c(:,end-filtpts+1:end);
cf(1,1:filtpts) = cf(1,1:filtpts) + 360.0;
cf(:,end-filtpts+1:end) = c(:,1:filtpts);
cf(1,end-filtpts+1:end) = cf(1,end-filtpts+1:end)-360.0;
cf(1,:) = filter_field(cf(1,:),filtpts,'-t');
cf(2,:) = filter_field(cf(2,:),filtpts,'-t');
c = cf(:,(filtpts+1):(filtpts+1+npts-1));
if (isobath == 100)
    inds = c(1,:)<-170;c(:,inds) = [];
    inds = c(1,:)>-78;c(:,inds) = [];
end

%Plot smoothed isobath:
hold on;
plot(c(1,:),c(2,:),'-','LineWidth',2,'color',[0 0.5 0]);

%Get x-y angle:
c(3,:) = zeros(size(c(1,:))); %x-y angle
cL = length(c(1,:));
for i = 2:(cL-1)
    dy = lat_to_km*(c(2,i+1)-c(2,i-1));
    dx = lat_to_km*(c(1,i+1)-c(1,i-1))*cos(pi*c(2,i)/180);
    c(3,i) = atan2(dy,dx)-pi;
end
c = c(:,2:(end-1));
c(3,c(3,:)<-pi) = c(3,c(3,:)<-pi)+2*pi;

cL = length(c(1,:));

%Angle = 0 => perpendicular to coast is north.
%Angle > 0 => perpendicular to coast is north-west.
%Angle < 0 => perpendicular to coast is north-east.
dx = zeros(length(c(1,:))-1,1);
for i=1:length(dx)
    dx(i) = Haversine(c(1,i),c(2,i),c(1,i+1),c(2,i+1));
end
dx = mean(dx)/1e3;

%Get segments:
sL = length(isp:sp:(cL-sp));

cseg = zeros(3,sL); %center (lon, lat, ang).
for i=1:sL
    cseg(1,i) = mean(c(1,(((i-1)*sp+1):(i*sp))+isp-1));
    cseg(2,i) = mean(c(2,(((i-1)*sp+1):(i*sp))+isp-1));
    cseg(3,i) = mean(c(3,(((i-1)*sp+1):(i*sp))+isp-1));
end

%Along-coast coordinate:
dc = zeros(sL,1);
for i=2:length(dc)
    dct = Haversine(cseg(1,i-1),cseg(2,i-1),cseg(1,i),cseg(2,i));
    dc(i) = dc(i-1)+dct;
end
dc = dc/1e3;

% $$$ hold on;
% $$$ plot(cseg(1,:),cseg(2,:),'+','LineWidth',1,'color',[0 0.5 0]);
% $$$ legend('Coast',sprintf('%03dm isobath',isobath),...
% $$$        sprintf('smoothed %03dm isobath',isobath),'Segment center positions');
legend('Coast',sprintf('%03dm isobath',isobath),...
       sprintf('smoothed %03dm isobath',isobath));

if (~plotting)
    close gcf;
end

end

