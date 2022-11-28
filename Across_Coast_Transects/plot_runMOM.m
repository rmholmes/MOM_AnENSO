% This function plots across coast transects from the MOM runs 

close all;
clear all;

Re = 6378000.0;
lat_to_km = 2*pi*Re/360.0/1e3;
zL = 50;

% MOM:
base = '/g/data/e14/rmh561/access-om2/archive/025deg_jra55_ryf/';

% 1/4deg:
grd_file = [base 'output019/ocean/ocean-2d-ht.nc'];
C_file = [base 'output000/ocean/ocean-3d-temp-1-monthly-mean-ym_1900_01.nc'];
D_file = [base 'output019/ocean/ocean-3d-temp-1-monthly-mean-ym_1919_01.nc'];
z = -ncread(C_file,'st_ocean');

mxind = 250;
[x_rho,y_rho] = ndgrid(ncread(grd_file,'xt_ocean'),ncread(grd_file,'yt_ocean'));
h = ncread(grd_file,'ht');
x_rho = x_rho(:,1:mxind);
y_rho = y_rho(:,1:mxind);
h = h(:,1:mxind);
mask = 1*isnan(h);
[xL,yL] = size(x_rho);
x_psi = avg(avg(x_rho,1),2);
y_psi = avg(avg(y_rho,1),2);


%Get smooth isobath and segments:
isobath = 1000;

plotting = 1;
plotnice = 0;
sp = 5; %Number of points to average angle and center point.
isp = 1; %Initial index to use
[c,dx,cseg,dc] = get_isobath_sections(x_rho,y_rho,mask,h,isobath,'MOM',plotting,sp,isp);

cL = length(c(1,:));
sL = length(cseg(1,:));

%Segment rotation options:
L = 2*sp*dx; %along-coast length of segment
W = 300; %off-coast distance of segment
Wm = 500; %on-coast distance of segment
Nw = round((W+Wm)/dx); %across-coast points to use
Nl = round(L/dx); %along-coast points to use

% Rotated segments:
for i=1:sL
    [lonr, latr, Corners,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,i));    
    if (~plotnice)
    hold on;
    plot([Corners(:,1); Corners(1,1)],[Corners(:,2); Corners(1,2)],'-m');    
    end
end
SECS = {};
names = {};

% 1000m isobath sections:
SECS{1} = [150:(150+22)] %hr
for i=SECS{1}
    [lonr, latr, Corners,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,i));
    if (~plotnice)
    hold on;
    plot([Corners(:,1); Corners(1,1)],[Corners(:,2); Corners(1,2)],'-g');
    end
end
names{1} = 'Bellingshausen'
if (plotnice)
    [lonr, latr, Corners1,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,SECS{1}(1)));
    [lonr, latr, Corners2,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,SECS{1}(end)));
    plot([Corners2(3,1) Corners1(4,1) Corners1(1,1) Corners2(2,1) Corners2(3,1)],...
         [Corners2(3,2) Corners1(4,2) Corners1(1,2) Corners2(2,2) Corners2(3,2)],'-g','linewidth',2);
    text(-130,-70,names{1},'color','g');
end


%% Single section plots:
    figure;
    set(gcf,'Position',[1         375        1918         584]);
    set(gcf,'defaulttextfontsize',15);
    set(gcf,'defaultaxesfontsize',15);
    poss = [0.0847    0.1134    0.3026    0.8150; ...
            0.5703    0.1134    0.3026    0.8150];
for tiii = 1:1
    clf;
    
    ti = [tiii tiii];

    files = {C_file,D_file}
    ziles = {z,z}

    for sss = 1:1
    seg = SECS{sss};

    titles = {['ACCESS-OM2 ' names{sss} ' Sector ' num2str(ti(1)+1968-1)]};

    [lonr, latr, Corners,ccd,lcd] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,seg(1)));    

    temp = zeros(length(ccd),zL,length(files),length(seg));
% $$$     dens = zeros(length(ccd),zL,length(files),length(seg));
    Z = zeros(length(ccd),zL,length(files),length(seg));
    Count = zeros(length(ccd),zL,length(files),length(seg));
    
    cnt = 0;
    for k=1:length(seg)
        [lonr, latr, Corners,ccd,lcd] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,seg(k)));    
        sprintf(['Doing segment %02d of ' names{sss}],k)
        for i = 1:length(files)
            
            if (i == 1)
                tiss = [1 1];
            else
                tiss = ti;
            end
            
            [T, tmp] = get_rotated_field(files{i},'temp','',tiss,lonr,latr,'t','3D',cseg(3,seg(k)));
            temp(:,:,i,k) = squeeze(nanmean(T,2))-273.15;
            T(~isnan(T)) = 1;T(isnan(T)) = 0;
            Count(:,:,i,k) = squeeze(sum(T,2));

            % If you have a vector field that you want to rotate:
% $$$             if (i ~= 2)
% $$$             [U, V] = get_rotated_field(files{i},'u','v',ti,lonr,latr,'c','3D',cseg(3,k));
% $$$             cc(:,:,i) = cc(:,:,i) + squeeze(nanmean(U,2));
% $$$             lc(:,:,i) = lc(:,:,i) + squeeze(nanmean(V,2));
% $$$             end
        end
        cnt = cnt+1;
    end
    temp = nanmean(temp,4);
% $$$     cc(:,:,:) = cc(:,:,:)/cnt;
% $$$     lc(:,:,:) = lc(:,:,:)/cnt;
% $$$     SSH(:,:,:) = SSH(:,:,:)/cnt;
    Count = sum(Count,4);
    Cfrac = Count/max(max(max(Count)));

    temp(Cfrac<0) = NaN;
    
    %Equalize NaNs:
    NaNs = repmat(isnan(lc(:,:,1)),[1 1 3]);
    temp(NaNs) = NaN;
    
    for i=1:length(ziles)
        Z(:,:,i) = repmat(ziles{i}',[length(ccd) 1]);
    end
    X = repmat(ccd',[1 zL]);

    %T control plots:
    %Settings:
    xlims = [-500 300];
    xtic  = 50;
    ylims = [-700 0];
    ytic = 500;
    cintR = [1020:0.1:1040];
% $$$     caxsT = [-1e10 -0.4:0.02:0.4 1e10];
    caxsT = [-1e10 -2:0.2:5 1e10];
    
    subplot(1,2,sss);
    
    contourf(X,Z(:,:,i),temp(:,:,2),caxsT,'linestyle','none');
    hold on;
    contour(X,Z(:,:,i),Cfrac(:,:,1),[0.9 0.9],'-k','LineWidth',2);
    contour(X,Z(:,:,i),temp(:,:,1),[-1e10 -2:0.25:2 1e10],'-','color',[0.5 0.5 0.5]);
    caxis([caxsT(2) caxsT(end-1)]);
    cb = colorbar;
    ylabel(cb,'Temperature Anomaly ($^\circ$C)');
    ylim(ylims);
    xlim(xlims);
    set(gca,'color','k');
    ylabel('Depth (m)');    
    xlabel('Across-Coast Distance (km)');
    title(titles{1});
    set(gca,'Position',poss(sss,:));
    end
end

