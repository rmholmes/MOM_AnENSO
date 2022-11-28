% This function plots across coast transects from the MOM runs and
% SOSE

close all;
clear all;

Re = 6378000.0;
lat_to_km = 2*pi*Re/360.0/1e3;
zL = 50;

% MOM:
base = '/srv/ccrc/data03/z3500785/MOM_AnENSO/';

% 1/4deg:
grd_file = [base 'mom025.ocean_grid.nc'];
C_file = [base 'mom025.566.ocean.nc'];
D_file = [base 'temp.cat01to25.diff.nc'];
z = -ncread(C_file,'st_ocean');
% $$$ W_file = [base 'ocean.201to210.cat.ncea.0to24.nc'];
% $$$ D_file = [base 'ocean.201to210.cat.diff.ncea.0to24.nc'];
% $$$ time = ncread(C_file,'time');
% $$$ t1 = 7.3003e4-0.5;
% $$$ time = time-t1;
% $$$ tL = length(time);
%Approximate Depth array:
% $$$ DZ = squeeze(nanmean(nanmean(ncread(C_file,'dzt',[1 1 1 1],[1440 150 50 1]),1),2));
% $$$ z = -cumsum(DZ);

% $$$ %1/10deg:
% $$$ grd_file01 = [base 'kds75.ocean_grid.nc'];
% $$$ C_file01 = [base 'kds75.ocean.tsrho.518to521.ncra.nc'];
% $$$ z01 = -ncread(C_file01,'st_ocean');
% $$$ z01 = z01(1:zL);
% $$$ 
% $$$ % ACCESS-OM2-025:
% $$$ grd_fileA025 = [base 'access-om2-025.ocean_grid.nc'];
% $$$ C_fileA025 = [base 'access-om2-025.072to075.ncra.nc'];
% $$$ zA025 = -ncread(C_fileA025,'st_ocean');
% $$$ 
% $$$ % ACCESS-OM2-01:
% $$$ grd_fileA01 = [base 'access-om2-01.ocean_grid.nc'];
% $$$ C_fileA01 = [base 'access-om2-01.428to431.ncra.nc'];
% $$$ zA01 = -ncread(C_fileA01,'st_ocean');
% $$$ zA01 = zA01(1:zL);

mxind = 250;
x_rho = ncread(grd_file,'geolon_t');
y_rho = ncread(grd_file,'geolat_t');
h = ncread(grd_file,'ht');
x_rho = x_rho(:,1:mxind);
y_rho = y_rho(:,1:mxind);
h = h(:,1:mxind);
mask = 1*isnan(h);
[xL,yL] = size(x_rho);
x_psi = avg(avg(x_rho,1),2);
y_psi = avg(avg(y_rho,1),2);


%Get smooth isobath and segments:
% $$$ isobath = 1000;
isobath = 100;

plotting = 1;
plotnice = 1;
sp = 5; %Number of points to average angle and center point.
%sp = 20; %lr
isp = 1; %Initial index to use
[c,dx,cseg,dc] = get_isobath_sections(x_rho,y_rho,mask,h,isobath,'MOM',plotting,sp,isp);

cL = length(c(1,:));
sL = length(cseg(1,:));

%Segment rotation options:
L = 2*sp*dx; %along-coast length of segment
W = 300; %off-coast distance of segment
Wm = 200; %on-coast distance of segment
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

% $$$ % 1000m isobath sections:
% $$$ SECS{1} = [285:(285+22)] %hr
% $$$ for i=SECS{1}
% $$$     [lonr, latr, Corners,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,i));
% $$$     if (~plotnice)
% $$$     hold on;
% $$$     plot([Corners(:,1); Corners(1,1)],[Corners(:,2); Corners(1,2)],'-k');
% $$$     end
% $$$ end
% $$$ names{1} = 'East Ross'
% $$$ if (plotnice)
% $$$     [lonr, latr, Corners1,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,SECS{1}(1)));
% $$$     [lonr, latr, Corners2,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,SECS{1}(end)));
% $$$     plot([Corners2(3,1) Corners1(4,1) Corners1(1,1) Corners2(2,1) Corners2(3,1)],...
% $$$          [Corners2(3,2) Corners1(4,2) Corners1(1,2) Corners2(2,2) Corners2(3,2)],'-g','linewidth',2);
% $$$     text(-130,-70,names{1},'color','g');
% $$$ end

% $$$ SECS{2} = [206:222] %hr
% $$$ for i=SECS{2}
% $$$     [lonr, latr, Corners,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,i));    
% $$$     if (plotting)
% $$$     hold on;
% $$$     plot([Corners(:,1); Corners(1,1)],[Corners(:,2); Corners(1,2)],'-k');
% $$$     end
% $$$ end
% $$$ names{2} = 'Western Peninsula (NATCC Fig. 3)'
% $$$ text(-100,-20,names{1});
% $$$ 
% 100m isobath sections:
SECS{1} = [1:16];%285:(285+22)] %hr
for i=SECS{1}
    [lonr, latr, Corners,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,i));
    if (~plotnice)
    hold on;
    plot([Corners(:,1); Corners(1,1)],[Corners(:,2); Corners(1,2)],'-k');
    end
end
names{1} = 'Bellingshausen'
if (plotnice)
    [lonr, latr, Corners1,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,SECS{1}(1)));
    [lonr, latr, Corners2,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,SECS{1}(end)));
    plot([Corners2(3,1) Corners1(4,1) Corners1(1,1) Corners2(2,1) Corners2(3,1)],...
         [Corners2(3,2) Corners1(4,2) Corners1(1,2) Corners2(2,2) Corners2(3,2)],'-m','linewidth',2);
    text(-130,-70,names{1},'color','m');
end

% 100m isobath sections:
SECS{2} = [47:61];%285:(285+22)] %hr
for i=SECS{2}
    [lonr, latr, Corners,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,i));
    if (~plotnice)
    hold on;
    plot([Corners(:,1); Corners(1,1)],[Corners(:,2); Corners(1,2)],'-k');
    end
end
names{2} = 'Amundsen'
if (plotnice)
    [lonr, latr, Corners1,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,SECS{2}(1)));
    [lonr, latr, Corners2,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,SECS{2}(end)));
    plot([Corners2(3,1) Corners1(4,1) Corners1(1,1) Corners2(2,1) Corners2(3,1)],...
         [Corners2(3,2) Corners1(4,2) Corners1(1,2) Corners2(2,2) Corners2(3,2)],'-b','linewidth',2);
    text(-130,-70,names{2},'color','b');
end


%% Single section plots:
% $$$     %Time average:
% $$$     dayi = 125; %start day (actual)
% $$$     [tmp tii] = min(abs(time-dayi));
% $$$     dayf = 175; %end day (actual)
% $$$     [tmp tif] = min(abs(time-(dayf-5)));
% $$$     ti = [tii tif]; 
    figure;
    set(gcf,'Position',[1         375        1918         584]);
    set(gcf,'defaulttextfontsize',15);
    set(gcf,'defaultaxesfontsize',15);
    poss = [0.0847    0.1134    0.3026    0.8150; ...
            0.5703    0.1134    0.3026    0.8150];
for tiii = 1:50
    clf;
    
    ti = [tiii tiii];

    %Get and plot variables along segments:
% $$$     files = {C_file,C_file01};
% $$$     ziles = {z,z01};
% $$$     files = {C_file,C_fileA025};
% $$$     ziles = {z,zA025};
% $$$     files = {C_file01,C_fileA01};
% $$$     ziles = {z01,zA01};
    files = {C_file,D_file}
    ziles = {z,z}

    for sss = 1:2
% $$$     sss = 1;
    seg = SECS{sss};

% $$$     titles = {['MOM025 ' names{sss} ' Sector'],['MOM01 ' names{sss} ' Sector']};
% $$$     titles = {['MOM025 ' names{sss} ' Sector'],['ACCESS-OM2-025 ' names{sss} ' Sector']};
% $$$     titles = {['MOM-SIS01 ' names{sss} ' Sector'],['ACCESS-OM2-01 ' names{sss} ' Sector']};
    titles = {['MOM025 ' names{sss} ' Sector ' num2str(ti(1)+1968-1)]};

    [lonr, latr, Corners,ccd,lcd] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,seg(1)));    

    temp = zeros(length(ccd),zL,length(files),length(seg));
    dens = zeros(length(ccd),zL,length(files),length(seg));
    Z = zeros(length(ccd),zL,length(files),length(seg));
    Count = zeros(length(ccd),zL,length(files),length(seg));
% $$$     cc = zeros(length(ccd),zL,length(files),length(seg));
% $$$     lc = zeros(length(ccd),zL,length(files),length(seg));
% $$$     SSH  = zeros(length(ccd),1,length(files),length(seg));
% $$$     dep = zeros(length(ccd),zL);
    
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
            temp(:,:,i,k) = squeeze(nanmean(T,2));
            T(~isnan(T)) = 1;T(isnan(T)) = 0;
            Count(:,:,i,k) = squeeze(sum(T,2));
            if (i == 1)
                [T, tmp] = get_rotated_field(files{i},'pot_rho_0','',tiss,lonr,latr,'t','3D',cseg(3,seg(k)));
                dens(:,:,i,k) = squeeze(nanmean(T,2));
            else
                dens(:,:,i,k) = zeros(size(Count(:,:,i,k)));
            end

% $$$             if (i ~= 2)
% $$$             [U, V] = get_rotated_field(files{i},'u','v',ti,lonr,latr,'c','3D',cseg(3,k));
% $$$             cc(:,:,i) = cc(:,:,i) + squeeze(nanmean(U,2));
% $$$             lc(:,:,i) = lc(:,:,i) + squeeze(nanmean(V,2));
% $$$             end
% $$$             
% $$$             [ssh, tmp] = get_rotated_field(files{i},'sea_level','',ti,lonr,latr,'t','2D',cseg(3,k));
% $$$             SSH(:,:,i) = SSH(:,:,i) + squeeze(nanmean(ssh,2));
        end
        cnt = cnt+1;
    end
    temp = nanmean(temp,4);
    dens = nanmean(dens,4);
% $$$     dep(:,:) = dep(:,:)/cnt;
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
    xlims = [-200 300];
    xtic  = 50;
    ylims = [-700 0];
    ytic = 500;
    cintR = [1020:0.1:1040];
    caxsT = [-1e10 -0.4:0.02:0.4 1e10];
    
    subplot(1,2,sss);
    
    contourf(X,Z(:,:,i),temp(:,:,2),caxsT,'linestyle','none');
    hold on;
    contour(X,Z(:,:,i),Cfrac(:,:,1),[0.9 0.9],'-k','LineWidth',2);
    contour(X,Z(:,:,i),dens(:,:,1),cintR,'-k');
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
colormap(redblue(40));
export_fig(sprintf('Tanoms_Movie_%04d.png',tiii+1968-1),'-painters');
end

% $$$     %T-anomaly plots:
% $$$     %Settings:
% $$$     xlims = [-200 300];
% $$$     xtic  = 50;
% $$$     ylims = [-1000 0];
% $$$     ytic = 500;
% $$$     cintR = [1020:0.1:1040];
% $$$     caxsT = [-1e10 -2:0.2:2 1e10];
% $$$     
% $$$     figure;
% $$$     clf;
% $$$     set(gcf,'Position',[1         375        1918         584]);
% $$$     set(gcf,'defaulttextfontsize',15);
% $$$     set(gcf,'defaultaxesfontsize',15);
% $$$     poss = [0.05 0.1 0.25 0.8; ...
% $$$             0.36 0.1 0.25 0.8; ...
% $$$             0.672 0.1 0.25 0.8]
% $$$ % $$$     for i=1:2
% $$$         subplot(1,3,i);
% $$$         contourf(X,Z(:,:,i),temp(:,:,i)+(1-i)*273.15,caxsT,'linestyle','none');
% $$$         hold on;
% $$$         hold on;
% $$$         contour(X,Z(:,:,i),Cfrac(:,:,i),[0.9 0.9],'-k','LineWidth',2);
% $$$         contour(X,Z(:,:,i),dens(:,:,i),cintR,'-k');
% $$$         caxis([caxsT(2) caxsT(end-1)]);
% $$$         ylim(ylims);
% $$$         xlim(xlims);
% $$$ % $$$         cb = colorbar;
% $$$ % $$$         ylabel(cb,'Control Temperature ($^\circ$C)');
% $$$         set(gca,'color','k');
% $$$         ylabel('Depth (m)');    
% $$$         xlabel('Across-Coast Distance (km)');
% $$$         title(titles{sss});
% $$$         set(gca,'Position',poss(i,:));
% $$$     end
% $$$ colormap(redblue(40));

% $$$     caxsU = [-1e10 -0.1:0.005:0.1 1e10]*100;
% $$$     caxsUdif = [-1e10 -0.05:0.0025:0.05 1e10]*100;
% $$$     %    caxsUdif = [-1e10 -0.02:0.002:0.02 1e10]*100;
% $$$     %caxsUdif = [-1e10 -0.01:0.001:0.01 1e10]*100;
% $$$ 
% $$$     cintR = [1020:0.04:1040];
% $$$     cintRdifP = [0.004:0.004:0.2];
% $$$     cintRdifM = [-0.2:0.004:-0.004];
% $$$ % $$$     cintRdifP = [0.002:0.002:0.1];
% $$$ % $$$     cintRdifM = [-0.1:0.002:-0.002];
% $$$ % $$$     cintRdifP = [0.001:0.001:0.05];
% $$$ % $$$     cintRdifM = [-0.05:0.001:-0.001];
% $$$ 
% $$$ % $$$     cintT = [-2:0.5:5];
% $$$ % $$$     cintTdifP = [0.1:0.1:0.5];
% $$$ % $$$     cintTdifM = [-0.5:0.1:-0.1];
% $$$ 
% $$$     SSHl = [-160 -140];
% $$$     SSHldif = [-5 0];
% $$$     
% $$$     %    labels = {'Control','WIND','DIFF'};    
% $$$     % Control T
% $$$     subplot(5,2,[3 5]);
% $$$     %    ax = axes('Position',posm(ii,:));
% $$$     contourf(X,Z,temp(:,:,1),caxsT,'linestyle','none');
% $$$     hold on;
% $$$     contour(X,Z,dens(:,:,1),cintR,'-k');
% $$$     contour(X,Z,dens(:,:,2),cintR,'--k');
% $$$     caxis([caxsT(2) caxsT(end-1)]);
% $$$     ylim(ylims);
% $$$     xlim(xlims);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Control Temperature ($^\circ$C)');
% $$$     set(gca,'Position',[0.07    0.4553    0.3385    0.2970]);
% $$$     set(gca,'color','k');
% $$$     ylabel('Depth (m)');
% $$$     text(128.9,1450,[names{sss} ' days ' num2str(dayi) ' to ' num2str(dayf)],'FontSize',20);
% $$$     text(txtx,txty,{'Control (-), Perturbation (- -)',sprintf('%.2fkgm$^{-3}$ Isopycnals',cintR(2)-cintR(1))},'FontSize',12,'color','w');
% $$$     
% $$$     subplot(5,2,[7 9]);
% $$$     %    ax = axes('Position',posm(ii,:));
% $$$     contourf(X,Z,lc(:,:,1)*100,caxsU,'linestyle','none');
% $$$ % $$$     hold on;
% $$$ % $$$     contour(X,Z,temp(:,:,1),cintT,'-k');
% $$$ % $$$     contour(X,Z,temp(:,:,2),cintT,'--k');
% $$$     caxis([caxsU(2) caxsU(end-1)]);
% $$$     ylim(ylims);
% $$$     xlim(xlims);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Control Along-coast Velocity  (cms$^{-1}$)');
% $$$     set(gca,'Position',[0.07    0.1100    0.3385    0.2970]);
% $$$     set(gca,'color','k');
% $$$     xlabel('Across-shelf Distance (km)');
% $$$     ylabel('Depth (m)');
% $$$ % $$$     text(txtx,txty,'Control');
% $$$     
% $$$     subplot(5,2,1);
% $$$     plot(X(:,1),SSH(:,1,1)*100,'-k');
% $$$     hold on;
% $$$     plot(X(:,1),SSH(:,1,2)*100,'--k');
% $$$     ylim(SSHl);
% $$$     xlim(xlims);
% $$$     set(gca,'Position',[0.07    0.8007    0.3385    0.1243]);
% $$$     ylabel('SSH (cm)');
% $$$     lg = legend('Control','Perturbed');
% $$$     set(lg,'Position',[0.2567    0.8809    0.1104    0.0514]);
% $$$     
% $$$     % Diff T
% $$$     subplot(5,2,[4 6]);
% $$$     %    ax = axes('Position',posm(ii,:));
% $$$     contourf(X,Z,temp(:,:,3),caxsTdif,'linestyle','none');
% $$$     hold on;
% $$$     contour(X,Z,dens(:,:,3),cintRdifP,'-k');
% $$$     contour(X,Z,dens(:,:,3),cintRdifM,'--k');
% $$$     caxis([caxsTdif(2) caxsTdif(end-1)]);
% $$$     ylim(ylims);
% $$$     xlim(xlims);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Temperature difference ($^\circ$C)');
% $$$     set(gca,'Position',[0.56    0.4553    0.3385    0.2970]);
% $$$     set(gca,'color','k');
% $$$     text(txtx,txty,{'Positive (-), Negative (- -)',['Density Difference ' ...
% $$$                      'Contours']},'FontSize',12,'color','w');
% $$$     %    text(txtx,txty,{'Temperature and','Density Difference'});
% $$$             
% $$$     subplot(5,2,[8 10]);
% $$$     %    ax = axes('Position',posm(ii,:));
% $$$     contourf(X,Z,lc(:,:,3)*100,caxsUdif,'linestyle','none');
% $$$ % $$$     hold on;
% $$$ % $$$     contour(X,Z,temp(:,:,3),cintTdifP,'-k');
% $$$ % $$$     contour(X,Z,temp(:,:,3),cintTdifM,'--k');
% $$$     caxis([caxsUdif(2) caxsUdif(end-1)]);
% $$$     ylim(ylims);
% $$$     xlim(xlims);
% $$$     cb = colorbar;
% $$$     ylabel(cb,'Along-coast velocity difference (cms$^{-1}$)');
% $$$     set(gca,'Position',[0.56    0.1100    0.3385    0.2970]);
% $$$     set(gca,'color','k');
% $$$     xlabel('Across-shelf Distance (km)');
% $$$     %    text(txtx,txty,{'Difference'Along-coast Vel and','Temperature Difference'});
% $$$ 
% $$$     subplot(5,2,2);
% $$$     plot(X(:,1),SSH(:,1,3)*100,'-k');
% $$$     ylim(SSHldif);
% $$$     xlim(xlims);
% $$$     set(gca,'Position',[0.56    0.8007    0.3385    0.1243]);
% $$$     ylabel('SSH Anomaly (cm)');
% $$$     set(lg,'Position',[0.5883    0.9156    0.1085    0.0280]);
% $$$ colormap(redblue(40));

% $$$ %% X-Y Plots of Temperature difference:
% $$$ 
% $$$ load('10thdegree/Tdiff_50to200m.mat');
% $$$ domain = [-240 0 -78 -58];
% $$$ [tmp mnln] = min(abs(lon-domain(1)));
% $$$ [tmp mxln] = min(abs(lon-domain(2)));
% $$$ [tmp mnlt] = min(abs(lat-domain(3)));
% $$$ v[tmp mxlt] = min(abs(lat-domain(4)));
% $$$ 
% $$$ [LON,LAT] = ndgrid(lon(mnln:mxln),lat(mnlt:mxlt));
% $$$ 
% $$$ figure;
% $$$ set(gcf,'defaulttextfontsize',10);
% $$$ set(gcf,'defaultaxesfontsize',10);
% $$$ set(gcf,'Position',[1          36        1920         970]);
% $$$ subplot(3,1,1);
% $$$ pcolPlot(LON,LAT,temp41(mnln:mxln,mnlt:mxlt));
% $$$ set(gca,'color','k');
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ title('Year 41');
% $$$ caxis([-0.5 0.5]);
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Temperature Difference 50-200m ($^\circ$C)');
% $$$ 
% $$$ subplot(3,1,2);
% $$$ pcolPlot(LON,LAT,temp44(mnln:mxln,mnlt:mxlt));
% $$$ set(gca,'color','k');
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ title('Year 44');
% $$$ caxis([-1.5 1.5]);
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Temperature Difference 50-200m ($^\circ$C)');
% $$$ 
% $$$ subplot(3,1,3);
% $$$ pcolPlot(LON,LAT,temp45(mnln:mxln,mnlt:mxlt));
% $$$ set(gca,'color','k');
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ title('Year 45');
% $$$ caxis([-1.5 1.5]);
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Temperature Difference 50-200m ($^\circ$C)');

%% SOSE T section: 
theta_file = [base 'sose_theta.ltmm.nc'];
salt_file = [base 'sose_salt.ltmm.nc'];
zSOSE = ncread(theta_file,'depth');
zL = 42;
% $$$ SOSE_lon = ncread(theta_file,'longitude')-360;
% $$$ SOSE_lat = ncread(theta_file,'latitude');
% $$$ SOSE_time = ncread(theta_file,'time');
% $$$ [tmp SOSEyL] = min(abs(SOSE_lat+58));
% $$$ SOSE_lat = SOSE_lat(1:SOSEyL);
% $$$ [SOSEx_rho,SOSEy_rho] = meshgrid(SOSE_lon,SOSE_lat);
% $$$ SOSExL = length(SOSE_lon);
% $$$ SOSEtL = length(SOSE_time);

%Get and plot variables along segments:
% $$$ sss = 1;
seg = SECS{sss};
[lonr, latr, Corners,ccd,lcd] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,seg(1)));    
    
i=1;
ti = [1 12];
SOSEtemp = zeros(length(ccd),zL,length(seg));
SOSEsalt = zeros(length(ccd),zL,length(seg));
SOSEdens = zeros(length(ccd),zL,length(seg));
SOSEcont = zeros(length(ccd),zL,length(seg));
    
cnt = 0;
for k=1:length(seg)
    sprintf(['Doing segment %02d of ' names{sss}],k)
    [T, tmp] = get_rotated_field(theta_file,'temperature','',ti,lonr,latr,'t','3D',cseg(3,seg(k)));
    SOSEtemp(:,:,k) = squeeze(nanmean(T,2));
    T(~isnan(T)) = 1;T(isnan(T)) = 0;
    SOSEcont(:,:,k) = squeeze(sum(T,2));
    [T, tmp] = get_rotated_field(salt_file,'salt','',ti,lonr,latr,'t','3D',cseg(3,seg(k)));
    SOSEsalt(:,:,k) = squeeze(nanmean(T,2));
    cnt = cnt+1;
end
SOSEdens = sw_dens0(SOSEsalt,SOSEtemp);
SOSEtemp = nanmean(SOSEtemp,3);
SOSEdens = nanmean(SOSEdens,3);

SOSEcont = sum(SOSEcont,4);
SOSECfrac = SOSEcont/max(max(max(SOSEcont)));

SOSEtemp(SOSECfrac<0) = NaN;

SOSEZ = repmat(zSOSE',[length(ccd) 1]);
SOSEX = repmat(ccd',[1 zL]);

subplot(1,3,3);
contourf(SOSEX,SOSEZ,SOSEtemp,caxsT,'linestyle','none');
hold on;
contour(SOSEX,SOSEZ(:,:,1),SOSECfrac(:,:,1),[0.9 0.9],'-k','LineWidth',2);
contour(SOSEX,SOSEZ(:,:,1),SOSEdens(:,:,1),cintR,'-k');
caxis([caxsT(2) caxsT(end-1)]);
ylim(ylims);
xlim(xlims);
cb = colorbar;
ylabel(cb,'Control Temperature ($^\circ$C)');
set(gca,'color','k');
ylabel('Depth (m)');    
xlabel('Across-Coast Distance (km)');
title(['SOSE ' names{sss} ' Sector']);
set(gca,'Position',poss(3,:));

% $$$ figure;
% $$$ subplot(3,2,1);
% $$$ contourf(avg(X,1),avg(Z,1),diff(temp(:,:,1),[],1)./diff(X(:,:,1), ...
% $$$                                                   [],1),[-100 ...
% $$$                     -0.05:0.005:0.05 100],'linestyle','none');
% $$$ ylim(ylims);
% $$$ xlim(xlims);
% $$$ % $$$ xlabel('Across-shelf Distance (km)');
% $$$ ylabel('Depth (m)');
% $$$ title('MOM dTdx ($^\circ$Ckm$^{-1}$');
% $$$ caxis([-0.05 0.05]);
% $$$ colorbar;
% $$$ set(gca,'color','k');
% $$$ 
% $$$ subplot(3,2,2);
% $$$ contourf(avg(SOSEX,1),avg(SOSEZ,1),diff(SOSEtemp(:,:,1),[],1)./diff(SOSEX(:,:,1), ...
% $$$                                                   [],1),[-100 ...
% $$$                     -0.05:0.005:0.05 100],'linestyle','none');
% $$$ ylim(ylims);
% $$$ xlim(xlims);
% $$$ % $$$ xlabel('Across-shelf Distance (km)');
% $$$ ylabel('Depth (m)');
% $$$ title('SOSE dTdx ($^\circ$Ckm$^{-1}$)');
% $$$ caxis([-0.05 0.05]);
% $$$ colorbar;
% $$$ set(gca,'color','k');
% $$$ 
% $$$ SOSEtempI = interp2(SOSEX',SOSEZ',SOSEtemp',X,Z,'linear');
% $$$ subplot(3,2,5);
% $$$ contourf(avg(X,1),avg(Z,1),abs(diff(temp(:,:,1),[],1))./abs(diff(SOSEtempI(:,:,1), ...
% $$$                                                   [],1)),[-100 0:0.5:5 ...
% $$$                    100],'linestyle','none');
% $$$ ylim(ylims);
% $$$ xlim(xlims);
% $$$ xlabel('Across-shelf Distance (km)');
% $$$ ylabel('Depth (m)');
% $$$ title('Ratio $|dTdx|_{MOM}/|dTdx|_{SOSE}$');
% $$$ caxis([0 5]);
% $$$ colorbar;
% $$$ set(gca,'color','k');
% $$$ 
% $$$ 
% $$$ subplot(3,2,3);
% $$$ contourf(avg(X,2),avg(Z,2),diff(temp(:,:,1),[],2)./diff(Z(:,:,1), ...
% $$$                                                   [],2),[-100 ...
% $$$                     -4e-3:1e-4:4e-3 100],'linestyle','none');
% $$$ ylim(ylims);
% $$$ xlim(xlims);
% $$$ % $$$ xlabel('Across-shelf Distance (km)');
% $$$ ylabel('Depth (m)');
% $$$ title('MOM dTdz ($^\circ$Cm$^{-1}$)');
% $$$ caxis([-0.004 0.004]);
% $$$ colorbar;
% $$$ set(gca,'color','k');
% $$$ 
% $$$ subplot(3,2,4);
% $$$ contourf(avg(SOSEX,2),avg(SOSEZ,2),diff(SOSEtemp(:,:,1),[],2)./diff(SOSEZ(:,:,1), ...
% $$$                                                   [],2),[-100 ...
% $$$                     -4e-3:1e-4:4e-3 100],'linestyle','none');
% $$$ ylim(ylims);
% $$$ xlim(xlims);
% $$$ % $$$ xlabel('Across-shelf Distance (km)');
% $$$ ylabel('Depth (m)');
% $$$ title('SOSE dTdz');
% $$$ caxis([-0.004 0.004]);
% $$$ colorbar;
% $$$ set(gca,'color','k');
% $$$ 
% $$$ SOSEtempI = interp2(SOSEX',SOSEZ',SOSEtemp',X,Z,'linear');
% $$$ subplot(3,2,6);
% $$$ contourf(avg(X,2),avg(Z,2),abs(diff(temp(:,:,1),[],2))./abs(diff(SOSEtempI(:,:,1), ...
% $$$                                                   [],2)),[-100 0:0.5:5 ...
% $$$                     100],'linestyle','none');
% $$$ ylim(ylims);
% $$$ xlim(xlims);
% $$$ xlabel('Across-shelf Distance (km)');
% $$$ ylabel('Depth (m)');
% $$$ title('Ratio $|dTdz|_{MOM}/|dTdz|_{SOSE}$');
% $$$ caxis([0 5]);
% $$$ colorbar;
% $$$ set(gca,'color','k');




% $$$ 
% $$$ %% Yearly-averaged diff plots:
% $$$     ti = [4]; %years to plot
% $$$     
% $$$     %Get and plot variables along segments:
% $$$     file = D_file;
% $$$     sss = 1;
% $$$     seg = SECS{sss};    
% $$$     [lonr, latr, Corners,ccd,lcd] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,seg(1)));    
% $$$     
% $$$     temp = zeros(length(ccd),zL,length(ti));
% $$$     dens = zeros(length(ccd),zL,length(ti));
% $$$     cc = zeros(length(ccd),zL,length(ti));
% $$$     lc = zeros(length(ccd),zL,length(ti));
% $$$     SSH  = zeros(length(ccd),1,length(ti));
% $$$     
% $$$     for i = 1:length(ti)
% $$$         cnt = 0;
% $$$         for k=1:length(seg)
% $$$             sprintf(['Doing segment %02d of ' names{sss}],k)
% $$$ 
% $$$             [T, tmp] = get_rotated_field(file,'temp','',ti(i),lonr,latr,'t','3D',cseg(3,seg(k)));
% $$$             temp(:,:,i) = temp(:,:,i) + squeeze(nanmean(T,2));
% $$$             [T, tmp] = get_rotated_field(file,'pot_rho_0','',ti(i),lonr,latr,'t','3D',cseg(3,seg(k)));
% $$$             dens(:,:,i) = dens(:,:,i) + squeeze(nanmean(T,2));
% $$$ 
% $$$             [U, V] = get_rotated_field(file,'u','v',ti(i),lonr,latr,'c','3D',cseg(3,seg(k)));
% $$$             cc(:,:,i) = cc(:,:,i) + squeeze(nanmean(U,2));
% $$$             lc(:,:,i) = lc(:,:,i) + squeeze(nanmean(V,2));
% $$$             
% $$$             [ssh, tmp] = get_rotated_field(file,'sea_level','',ti(i),lonr,latr,'t','2D',cseg(3,seg(k)));
% $$$             SSH(:,:,i) = SSH(:,:,i) + squeeze(nanmean(ssh,2));
% $$$             
% $$$             cnt = cnt+1;
% $$$         end
% $$$     temp(:,:,i) = temp(:,:,i)/cnt;
% $$$     dens(:,:,i) = dens(:,:,i)/cnt;
% $$$     cc(:,:,i) = cc(:,:,i)/cnt;
% $$$     lc(:,:,i) = lc(:,:,i)/cnt;
% $$$     SSH(:,:,i) = SSH(:,:,i)/cnt;
% $$$     end
% $$$ 
% $$$     Z = repmat(z',[length(ccd) 1]);
% $$$     X = repmat(ccd',[1 zL]);
% $$$ 
% $$$     %Equalize NaNs:
% $$$     NaNs = repmat(isnan(lc(:,:,1)),[1 1 3]);
% $$$     temp(NaNs) = NaN;
% $$$     
% $$$     %Settings:
% $$$     txtx = -195;
% $$$     txty = -1844;
% $$$     xlims = [-200 210];
% $$$     xtic  = 100;
% $$$     ylims = [-2000 0];
% $$$     ytic = 500;
% $$$     
% $$$     caxsTdif = [-1e10 -3:0.05:3 1e10];
% $$$     caxsUdif = [-1e10 -0.075:0.0025:0.075 1e10]*100;
% $$$ 
% $$$     cintRdifP = [0.008:0.008:0.8];
% $$$     cintRdifM = [-0.8:0.008:-0.008];
% $$$     
% $$$     SSHldif = [-7 0];
% $$$     
% $$$     %figure;
% $$$     clf;
% $$$     %    set(gcf,'Position',[1          36        1920         970]);
% $$$     set(gcf,'Position',[3 40 1904 963]);
% $$$     set(gcf,'defaulttextfontsize',10);
% $$$     set(gcf,'defaultaxesfontsize',10);
% $$$     
% $$$     for i = 1:length(ti)
% $$$ 
% $$$     subplot(5,length(ti),length(ti)*[1 2]+i);
% $$$     contourf(X,Z,temp(:,:,i),caxsTdif,'linestyle','none');
% $$$     hold on;
% $$$     contour(X,Z,dens(:,:,i),cintRdifP,'-k');
% $$$     contour(X,Z,dens(:,:,i),cintRdifM,'--k');
% $$$     caxis([caxsTdif(2) caxsTdif(end-1)]);
% $$$     ylim(ylims);
% $$$     xlim(xlims);
% $$$     set(gca,'color','k');
% $$$     text(txtx,txty,{'Positive (-), Negative (- -)',['Density Difference ' ...
% $$$                      'Contours']},'FontSize',10,'color','w');
% $$$     if (i==length(ti))
% $$$         pos = get(gca,'Position');
% $$$         cb = colorbar;
% $$$         ylabel(cb,'Temperature difference ($^\circ$C)');
% $$$         set(gca,'Position',pos);
% $$$     end
% $$$ 
% $$$     subplot(5,length(ti),length(ti)*[3 4]+i);
% $$$     contourf(X,Z,lc(:,:,i)*100,caxsUdif,'linestyle','none');
% $$$     caxis([caxsUdif(2) caxsUdif(end-1)]);
% $$$     ylim(ylims);
% $$$     xlim(xlims);
% $$$     set(gca,'color','k');
% $$$     xlabel('Distance from 1000m isobath (km)');
% $$$     if (i==length(ti))
% $$$         pos = get(gca,'Position');
% $$$         cb = colorbar;
% $$$         ylabel(cb,'Along-coast velocity difference (cms$^{-1}$)');
% $$$         set(gca,'Position',pos);
% $$$     end
% $$$ 
% $$$     subplot(5,length(ti),i);
% $$$     plot(X(:,1),SSH(:,1,i)*100,'-k');
% $$$     ylim(SSHldif);
% $$$     xlim(xlims);
% $$$     ylabel('SSH anomaly (m)');
% $$$     title(['Year ' num2str(i)],'FontSize',15);
% $$$     end
% $$$     colormap(redblue(40));
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ %% SSH and SSH gradient: 
% $$$ dat_file = D_file;
% $$$ 
% $$$ dayi = 5; %start day (actual)
% $$$ [tmp tii] = min(abs(time-dayi));
% $$$ dayf = 15; %end day (actual)
% $$$ [tmp tif] = min(abs(time-(dayf-5)));
% $$$ dstr = [num2str(dayi) '-' num2str(dayf)];
% $$$ 
% $$$ [xL,yL] = size(x_rho);
% $$$ 
% $$$ zeta = ncread(dat_file,'sea_level',[1 1 tii],[xL yL tif-tii+1]);
% $$$ %zeta = ncread(dat_file,'zeta',[1 1 tii],[xL yL tif-tii+1]);
% $$$ zeta = mean(zeta,3);
% $$$ 
% $$$ %Smoothing:
% $$$ zeta = filter_field(zeta,5,'-s');
% $$$ 
% $$$ %BT velocity:
% $$$ u = mean(nanmean(ncread(dat_file,'u',[1 1 1 tii],[xL yL zL tif-tii+1]),3),4);
% $$$ v = mean(nanmean(ncread(dat_file,'v',[1 1 1 tii],[xL yL zL tif-tii+1]),3),4);
% $$$ ubar = sqrt(u.^2+v.^2);
% $$$ ubar = filter_field(ubar,3,'-s');
% $$$ 
% $$$ % $$$ figure;
% $$$ %set(gcf,'Position',[15         394        1838         556]);
% $$$ set(gcf,'Position',[30 40 1838 963]);
% $$$ clf;
% $$$ 
% $$$ % $$$ SSHcont = [-0.2:0.001:0.1];
% $$$ % $$$ gSSHcont = [-1e20 0:0.02e-4:0.2e-4 1e20];
% $$$ % $$$ gSSHaxs = [0 0.2e-4];
% $$$ 
% $$$ SSHcont = [-0.2:0.001:0.1];
% $$$ gSSHcont = [-1e20 0:0.05e-4:0.5e-4 1e20];
% $$$ gSSHaxs = [0 0.5e-4];
% $$$ 
% $$$ BTvelcont = [-1e20 0:0.05:5 1e20];
% $$$ BTvelaxs = [0 0.5];
% $$$ 
% $$$ %SSH gradient:
% $$$ x = lat_to_km*cos(pi*y_rho/180).*x_rho;
% $$$ y = lat_to_km*y_rho;
% $$$ gradSSH = sqrt(avg(diff(zeta,[],1)./diff(x,[],1),2).^2+...
% $$$           avg(diff(zeta,[],2)./diff(y,[],2),1).^2);
% $$$ 
% $$$ %Contour figure setup:
% $$$ set(gcf,'Position',[301          40        1102         963]);
% $$$ set(gcf,'defaulttextfontsize',15);
% $$$ set(gcf,'defaultaxesfontsize',15);
% $$$ subplot(6,2,1);
% $$$ set(gca,'Position',[0.0699    0.7601    0.8657    0.2118]);
% $$$ % $$$ %gSSH:
% $$$ % $$$ gSSH = zeros(size(x_rho));
% $$$ % $$$ gSSH(2:(end-1),2:(end-1)) = avg(avg(gradSSH,1),2);
% $$$ % $$$ gSSH(isnan(gSSH)) = 0.0;
% $$$ % $$$ gSSH(mask == 1) = NaN;
% $$$ % $$$ contourf(x_rho,y_rho,gSSH,gSSHcont,'linestyle','none');
% $$$ %BT velocity:
% $$$ contourf(x_rho,y_rho,ubar*100,BTvelcont,'linestyle','none');
% $$$ hold on;
% $$$ contour(x_rho,y_rho,zeta,SSHcont,'-','color',[0.6 0.6 0.6]);
% $$$ xlim([-240 0]);
% $$$ ylim([-78 -55]);
% $$$ %caxis(gSSHaxs);
% $$$ caxis(BTvelaxs);
% $$$ cmap = colormap(redblue(40));
% $$$ cmap = cmap(21:end,:);
% $$$ colormap(cmap);
% $$$ 
% $$$ [c,tmp,tmp,tmp] = get_isobath_sections(x_rho,y_rho,mask,h,2500,'MOM',0,sp,isp);
% $$$ plot(c(1,:),c(2,:),'-b');
% $$$ [c,dx,cseg,dc] = get_isobath_sections(x_rho,y_rho,mask,h,1000,'MOM',0,sp,isp);
% $$$ plot(c(1,:),c(2,:),'-b','Linewidth',1);
% $$$ xlabel('Longitude ($^\circ$E)');
% $$$ ylabel('Latitude ($^\circ$N)');
% $$$ title(['MOM SSH anomaly (0.1cm contours) and gradient (color), ' ...
% $$$        'days ' dstr]);
% $$$ colorbar;
% $$$ 
% $$$ txpos = [-97.0793  -75.2795; ...
% $$$          -238.2584  -70.9559; ...
% $$$          -134.7081  -77.2155; ...
% $$$          -25.7465  -77.0983];
% $$$ %for kk=1:length(SECS)
% $$$ %for kk=2:4
% $$$ for kk=8:8
% $$$ out = zeros(1,2);
% $$$ in = zeros(1,2);
% $$$ cnt = 1;
% $$$ for i=SECS{kk}
% $$$     [lonr, latr, Corners,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,i));    
% $$$     out(cnt,1) = Corners(1,1);
% $$$     out(cnt,2) = Corners(1,2);
% $$$     in(cnt,1) = Corners(3,1);
% $$$     in(cnt,2) = Corners(3,2);
% $$$     cnt = cnt+1;
% $$$ end
% $$$ 
% $$$ lx = [out(:,1); flipud(in(:,1)); out(1,1)];
% $$$ ly = [out(:,2); flipud(in(:,2)); out(1,2)];
% $$$ hold on;
% $$$ plot(lx,ly,'-g','LineWidth',2);
% $$$ %text(txpos(kk,1),txpos(kk,2),names{kk},'color','g','FontSize',12);
% $$$ end
% $$$ set(gca,'color','k');
% $$$ set(gca,'Position',[0.0699    0.7165    0.84    0.2118]);
% $$$ text(-239.7775,-56.6471,'(a)','BackgroundColor','w','margin',0.1,'FontSize',15);
% $$$ 
% $$$ %Put on number labels for cross-coast plots:
% $$$ dist = [0:1000:10000];
% $$$ [tmp x2] = min(abs(cseg(1,:))-0);
% $$$ hand = zeros(size(dist));
% $$$ for i=1:length(dist)
% $$$     [tmp ind] = min(abs(dc-dc(x2)-dist(i)));
% $$$     %    plot(cseg(1,ind),cseg(2,ind),'+','color',[0 0.5 0]);
% $$$     hand(i) = text(cseg(1,ind),cseg(2,ind)+1.5,num2str(dist(i)/1000), ...
% $$$                    'color','k','HorizontalAlignment','center','verticalAlignment','middle');
% $$$ end
% $$$ 
% $$$ lcdim = dc-dc(x2);
% $$$ 
% $$$ hand = copyobj(gca,gcf);
% $$$ axes(hand);
% $$$ colorbar;
% $$$ set(gca,'Position',[0.0699 0.4 0.84 0.2118]);
% $$$ out = findobj(gca);
% $$$ delete(out(end));
% $$$ delete(out(end-1));
% $$$ 
% $$$ %ROMS:
% $$$ ROMSgrd_file = [base 'ocean_grd025.nc'];
% $$$ ROMSdat_file = [base 'ocean_his.nc'];
% $$$ ROMStime = ncread(ROMSdat_file,'ocean_time')/86400;
% $$$ 
% $$$ ROMSx_rho = ncread(ROMSgrd_file,'lon_rho');
% $$$ ROMSy_rho = ncread(ROMSgrd_file,'lat_rho');
% $$$ ROMSmask = ncread(ROMSgrd_file,'mask_rho');
% $$$ [xL,yL] = size(ROMSx_rho);
% $$$ 
% $$$ dayi = 5;
% $$$ dayf = 15;
% $$$ [tmp tii] = min(abs(ROMStime-dayi));
% $$$ [tmp tif] = min(abs(ROMStime-(dayf-5)));
% $$$ 
% $$$ zeta = ncread(ROMSdat_file,'zeta',[1 1 tii],[xL yL tif-tii+1]);
% $$$ zeta = mean(zeta,3);
% $$$ zeta = filter_field(zeta,5,'-s');
% $$$ 
% $$$ ubar = ncread(ROMSdat_file,'ubar',[1 1 tii],[xL-1 yL tif-tii+1]);
% $$$ ubar = mean(ubar,3);
% $$$ vbar = ncread(ROMSdat_file,'vbar',[1 1 tii],[xL yL-1 tif-tii+1]);
% $$$ vbar = mean(vbar,3);
% $$$ uvbar = zeros(size(zeta));
% $$$ uvbar(2:(end-1),2:(end-1)) = sqrt(avg(ubar(:,2:(end-1)),1).^2+avg(vbar(2:(end-1),:),2).^2);
% $$$ uvbar = filter_field(uvbar,3,'-s');
% $$$ 
% $$$ 
% $$$ % $$$ x = lat_to_km*cos(pi*ROMSy_rho/180).*ROMSx_rho;
% $$$ % $$$ y = lat_to_km*ROMSy_rho;
% $$$ % $$$ gradSSH = sqrt(avg(diff(zeta,[],1)./diff(x,[],1),2).^2+...
% $$$ % $$$           avg(diff(zeta,[],2)./diff(y,[],2),1).^2);
% $$$ % $$$ 
% $$$ % $$$ gSSH = zeros(size(ROMSx_rho));
% $$$ % $$$ gSSH(2:(end-1),2:(end-1)) = avg(avg(gradSSH,1),2);
% $$$ % $$$ gSSH(isnan(gSSH)) = 0.0;
% $$$ % $$$ gSSH(ROMSmask == 0) = NaN;
% $$$ %h1 = contourf(ROMSx_rho-360,ROMSy_rho,gSSH,gSSHcont,'linestyle','none');
% $$$ h1 = contourf(ROMSx_rho-360,ROMSy_rho,uvbar*100,BTvelcont,'linestyle','none');
% $$$ hold on;
% $$$ h2 = contour(ROMSx_rho-360,ROMSy_rho,zeta,SSHcont,'-','color',[0.6 0.6 0.6]);
% $$$ out = findobj(gca);
% $$$ uistack(out(2),'bottom');
% $$$ uistack(out(3),'bottom');
% $$$ xlim([-240 0]);
% $$$ ylim([-78 -55]);
% $$$ %caxis(gSSHaxs);
% $$$ caxis(BTvelaxs);
% $$$ title(['Barotropic model SSH anomaly (0.1cm contours) and gradient ' ...
% $$$        '(color), days ' dstr]);
% $$$ uistack(out(3),'bottom');
% $$$ text(-239.7775,-56.6471,'(b)','BackgroundColor','w','margin',0.1,'FontSize',15);
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ %%% Along Isobath plots:
% $$$ 
% $$$ depthi = 200;
% $$$ [tmp dii] = min(abs(z + depthi));
% $$$ %depthi = round(-z(dii));
% $$$ 
% $$$ depthf = 700;
% $$$ [tmp dif] = min(abs(z + depthf));
% $$$ %depthf = round(-z(dif));
% $$$ 
% $$$ %Control Temperature (total from first year):
% $$$ C_file = [base 'cntrl.201to210.cat.ncea.0to24.nc'];
% $$$ temp1 = sum(ncread(C_file,'temp',[1 1 dii 1],[xL yL dif-dii+1 25]),4);
% $$$ C_file = [base 'cntrl.201to210.cat.ncea.25to49.nc'];
% $$$ temp2 = sum(ncread(C_file,'temp',[1 1 dii 1],[xL yL dif-dii+1 25]),4);
% $$$ C_file = [base 'cntrl.201to210.cat.ncea.50to72.nc'];
% $$$ temp3 = sum(ncread(C_file,'temp',[1 1 dii 1],[xL yL dif-dii+1 23]),4);
% $$$ temp = nanmean((temp1+temp2+temp3)/(25+25+23),3);
% $$$ clear temp1 temp2 temp3;
% $$$ 
% $$$ %MOM SSH:
% $$$ D_file = [base 'ocean.201to210.cat.diff.ncea.0to24.nc'];
% $$$ dayi = 60; %start day (actual)
% $$$ dayf = 120; %end day (actual)
% $$$ time = ncread(D_file,'time')-7.3003e4+0.5;
% $$$ [tmp tii] = min(abs(time-dayi));
% $$$ [tmp tif] = min(abs(time-(dayf-5)));
% $$$ zeta = ncread(D_file,'sea_level',[1 1 tii],[xL yL tif-tii+1]);
% $$$ zeta = mean(zeta,3);
% $$$ 
% $$$ %MOM BT velocity:
% $$$ MOMubar = mean(nanmean(ncread(D_file,'u',[1 1 1 tii],[xL yL zL tif-tii+1]),3),4);
% $$$ MOMvbar = mean(nanmean(ncread(D_file,'v',[1 1 1 tii],[xL yL zL tif-tii+1]),3),4);
% $$$ x_c = ncread(D_file,'geolon_c');
% $$$ y_c = ncread(D_file,'geolat_c');
% $$$ 
% $$$ %ROMS SSH:
% $$$ dayi = 5;
% $$$ dayf = 15;
% $$$ ROMSgrd_file = [base 'ocean_grd025.nc'];
% $$$ ROMSdat_file = [base 'ocean_his.nc'];
% $$$ ROMStime = ncread(ROMSdat_file,'ocean_time')/86400;
% $$$ ROMSx_rho = ncread(ROMSgrd_file,'lon_rho')-360;
% $$$ ROMSy_rho = ncread(ROMSgrd_file,'lat_rho');
% $$$ ROMSmask = ncread(ROMSgrd_file,'mask_rho');
% $$$ [ROMSxL,ROMSyL] = size(ROMSx_rho);
% $$$ [tmp tii] = min(abs(ROMStime-dayi));
% $$$ [tmp tif] = min(abs(ROMStime-(dayf-5)));
% $$$ ROMSzeta = ncread(ROMSdat_file,'zeta',[1 1 tii],[ROMSxL ROMSyL tif-tii+1]);
% $$$ ROMSzeta = mean(ROMSzeta,3);
% $$$ 
% $$$ %ROMS BT velocity:
% $$$ ROMSubar = ncread(ROMSdat_file,'ubar',[1 1 tii],[ROMSxL-1 ROMSyL tif-tii+1]);
% $$$ ROMSubar = mean(ROMSubar,3);
% $$$ ROMSvbar = ncread(ROMSdat_file,'vbar',[1 1 tii],[ROMSxL ROMSyL-1 tif-tii+1]);
% $$$ ROMSvbar = mean(ROMSvbar,3);
% $$$ ROMSubar = avg(ROMSubar,2);
% $$$ ROMSvbar = avg(ROMSvbar,1);
% $$$ ROMSx_psi = avg(avg(ROMSx_rho,1),2);
% $$$ ROMSy_psi = avg(avg(ROMSy_rho,1),2);
% $$$ 
% $$$ %MOM Temperature Anomaly (total from first year):
% $$$ D_file = [base 'ocean.201to210.cat.diff.ncea.0to24.nc'];
% $$$ temp1 = sum(ncread(D_file,'temp',[1 1 dii 1],[xL yL dif-dii+1 25]),4);
% $$$ D_file = [base 'ocean.201to210.cat.diff.ncea.25to49.nc'];
% $$$ temp2 = sum(ncread(D_file,'temp',[1 1 dii 1],[xL yL dif-dii+1 25]),4);
% $$$ D_file = [base 'ocean.201to210.cat.diff.ncea.50to72.nc'];
% $$$ temp3 = sum(ncread(D_file,'temp',[1 1 dii 1],[xL yL dif-dii+1 23]),4);
% $$$ tempD = nanmean((temp1+temp2+temp3)/(25+25+23),3);
% $$$ % $$$ clear temp1 temp2 temp3;
% $$$ % $$$ D_file = [base 'ocean.201to210.cat.diff.ncea.50to72.nc'];
% $$$ % $$$ dayi = 300; %start day (actual)
% $$$ % $$$ dayf = 360; %end day (actual)
% $$$ % $$$ time = ncread(D_file,'time')-7.3003e4+0.5;
% $$$ % $$$ [tmp tii] = min(abs(time-dayi));
% $$$ % $$$ [tmp tif] = min(abs(time-(dayf-5)));
% $$$ % $$$ tempD = ncread(D_file,'temp',[1 1 dii tii],[xL yL dif-dii+1 tif-tii+1]);
% $$$ % $$$ tempD = mean(nanmean(tempD,3),4);
% $$$ 
% $$$ %N2:
% $$$ % $$$ rho = ncread(C_file,'pot_rho_0',[1 1 dii-1 tii],[xL yL dif-dii+3 tif-tii+1]);
% $$$ % $$$ rho = mean(rho,4);
% $$$ % $$$ dzt = ncread(C_file,'dzt',[1 1 dii-1 tii],[xL yL dif-dii+3 tif-tii+1]);
% $$$ % $$$ dzt = mean(dzt,4);
% $$$ % $$$ N2 = mean(avg(9.81/1025*diff(rho,[],3)./avg(dzt,3),3),3);
% $$$ 
% $$$ %Bathymetry:
% $$$ h = ncread(grd_file,'ht',[1 1],[xL yL]);
% $$$ 
% $$$ %SOSE Temperature:
% $$$ theta_file = [base 'sose_theta.ltmm.nc'];
% $$$ salt_file = [base 'sose_salt.ltmm.nc'];
% $$$ SOSE_lon = ncread(theta_file,'longitude')-360;
% $$$ SOSE_lat = ncread(theta_file,'latitude');
% $$$ SOSE_depth = ncread(theta_file,'depth');
% $$$ SOSE_time = ncread(theta_file,'time');
% $$$ [tmp SOSEyL] = min(abs(SOSE_lat+58));
% $$$ SOSE_lat = SOSE_lat(1:SOSEyL);
% $$$ [SOSEx_rho,SOSEy_rho] = meshgrid(SOSE_lon,SOSE_lat);
% $$$ SOSExL = length(SOSE_lon);
% $$$ SOSEtL = length(SOSE_time);
% $$$ [tmp dii] = min(abs(SOSE_depth + depthi));
% $$$ [tmp dif] = min(abs(SOSE_depth + depthf));
% $$$ SOSE_theta = mean(nanmean(ncread(theta_file,'temperature',[1 1 dii 1],[SOSExL SOSEyL dif-dii+1 SOSEtL]),3),4);
% $$$ 
% $$$ %Get lon-lat rotated:
% $$$ SSHr = zeros(Nw+1,sL);
% $$$ hr = zeros(Nw+1,sL);
% $$$ Tr = zeros(Nw+1,sL);
% $$$ TDr = zeros(Nw+1,sL);
% $$$ BTr = zeros(Nw+1,sL);
% $$$ ROMSSSHr = zeros(Nw+1,sL);
% $$$ ROMSBTr = zeros(Nw+1,sL);
% $$$ TSOSEr = zeros(Nw+1,sL);
% $$$ 
% $$$ % $$$ N2r = zeros(Nw+1,sL);
% $$$ 
% $$$ W = 400;
% $$$ for k=1:sL
% $$$     k
% $$$     [lonr, latr, Corners,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,k));    
% $$$     
% $$$ % $$$     ROMSu = interp2(ROMSx_psi',ROMSy_psi',ROMSubar',lonr,latr,'linear');
% $$$ % $$$     ROMSu = nanmean(ROMSu,2);
% $$$ % $$$     ROMSv = interp2(ROMSx_psi',ROMSy_psi',ROMSvbar',lonr,latr,'linear');
% $$$ % $$$     ROMSv = nanmean(ROMSv,2);
% $$$ % $$$     ROMSBTr(:,k) = cos(atan2(ROMSv,ROMSu)-cseg(3,k)).*sqrt(ROMSu.^2+ROMSv.^2);
% $$$ 
% $$$     MOMu = interp2(x_c',y_c',MOMubar',lonr,latr,'linear');
% $$$     MOMu = nanmean(MOMu,2);
% $$$     MOMv = interp2(x_c',y_c',MOMvbar',lonr,latr,'linear');
% $$$     MOMv = nanmean(MOMv,2);
% $$$     BTr(:,k) = cos(atan2(MOMv,MOMu)-cseg(3,k)).*sqrt(MOMu.^2+MOMv.^2);
% $$$ 
% $$$     SSHrot = interp2(x_rho',y_rho',zeta',lonr,latr,'linear');
% $$$     SSHr(:,k) = nanmean(SSHrot,2);
% $$$ % $$$     ROMSSSHrot = interp2(ROMSx_rho',ROMSy_rho',ROMSzeta',lonr,latr,'linear');
% $$$ % $$$     ROMSSSHr(:,k) = nanmean(ROMSSSHrot,2);
% $$$ % $$$     hrot = interp2(x_rho',y_rho',h',lonr,latr,'linear');
% $$$ % $$$     hr(:,k) = nanmean(hrot,2);
% $$$ % $$$     Trot = interp2(x_rho',y_rho',temp',lonr,latr,'linear');
% $$$ % $$$     Tr(:,k) = nanmean(Trot,2);
% $$$ % $$$     TDrot = interp2(x_rho',y_rho',tempD',lonr,latr,'linear');
% $$$ % $$$     TDr(:,k) = nanmean(TDrot,2);
% $$$ % $$$     TSOSErot = interp2(SOSEx_rho,SOSEy_rho,SOSE_theta',lonr,latr,'linear');
% $$$ % $$$     TSOSEr(:,k) = nanmean(TSOSErot,2);
% $$$ 
% $$$ % $$$     N2rot = interp2(x_rho',y_rho',N2',lonr,latr,'linear');
% $$$ % $$$     N2r(:,k) = nanmean(N2rot,2);
% $$$ end
% $$$ 
% $$$ figure;
% $$$ clf;
% $$$ set(gcf,'Position',[354          28        1035         975]);
% $$$ set(gcf,'defaulttextfontsize',10);
% $$$ set(gcf,'defaultaxesfontsize',10);
% $$$ set(gcf,'color','w');
% $$$ 
% $$$ poss = zeros(6,4);
% $$$ poss(1,:) = [0.1034    0.05    0.7372    0.1254];
% $$$ for i=2:6
% $$$     poss(i,:) = poss(i-1,:);
% $$$     poss(i,2) = poss(i-1,2)+0.16;
% $$$ end
% $$$ 
% $$$ %xpts:
% $$$ [tmp x1] = min(abs(cseg(1,:)+240));
% $$$ [tmp x2] = min(abs(cseg(1,:))-0);
% $$$ 
% $$$ lcdim = dc-dc(x2);
% $$$ 
% $$$ W = 210;
% $$$ [X,Y] = meshgrid(lcdim,cc);
% $$$ 
% $$$ secpl = [6];
% $$$ 
% $$$ %Regions:
% $$$ xreg = zeros(5,length(secpl));
% $$$ yreg = zeros(5,length(secpl));
% $$$ for i=1:length(secpl)
% $$$     xreg(:,i) = [X(1,SECS{secpl(i)}(1)); ...
% $$$                  X(1,SECS{secpl(i)}(end)); ...
% $$$                  X(1,SECS{secpl(i)}(end)); ...
% $$$                  X(1,SECS{secpl(i)}(1)); ...
% $$$                  X(1,SECS{secpl(i)}(1))];
% $$$     yreg(:,i) = [-Wm; -Wm; W; W; -Wm];
% $$$ end
% $$$ 
% $$$ %SSH Anomaly average:
% $$$ SSHmom = nanmean(SSHr(:,SECS{6}(1):SECS{6}(end)),2);
% $$$ SSHroms = nanmean(ROMSSSHr(:,SECS{6}(1):SECS{6}(end)),2);
% $$$ BTroms = nanmean(ROMSBTr(:,SECS{6}(1):SECS{6}(end)),2);
% $$$ BT = nanmean(BTr(:,SECS{6}(1):SECS{6}(end)),2);
% $$$ 
% $$$ subplot(6,1,1);
% $$$ contourf(X,Y,hr,[0:500:4000 10000],'-k','linestyle','none');
% $$$ hold on;
% $$$ contour(X,Y,hr,[1000 1000],'-k');
% $$$ contour(X,Y,hr,[2500 2500],'-k');
% $$$ plot(xreg,yreg,'-g','LineWidth',2);
% $$$ ylabel({'   Across-shelf ','Distance (km)'});
% $$$ cb = colorbar;
% $$$ %ylabel(cb,'Bathymetry (m)');
% $$$ caxis([0 4000]);
% $$$ set(gca,'xdir','reverse');
% $$$ xlim([X(1,x2) X(1,x1)]);
% $$$ set(gca,'color','k');
% $$$ txpos = [2149.6 347.7; ...
% $$$          978.3 245.5; ...
% $$$          1714.4 352.3; ...
% $$$          3020.1 347.7];
% $$$ % $$$ for i=1:length(secpl)
% $$$ % $$$     text(2300,250,names{i},'color','g','FontSize',12);
% $$$ % $$$ end
% $$$ %text(X(1,x1)+35,350,'(a)','BackgroundColor','w','FontSize',12,'Margin',0.1);
% $$$ set(gca,'Position',poss(end,:));
% $$$ set(gca,'xticklabel',[]);
% $$$ 
% $$$ subplot(6,1,2);
% $$$ contourf(X,Y,BTr*100,[-1e20 -1:0.025:1 1e20],'-k','linestyle','none');
% $$$ contourf(X,Y,BTr*100,[-1e20 -5:0.1:5 1e20],'-k','linestyle','none');
% $$$ hold on;
% $$$ contour(X,Y,SSHr*100,[-1e20 -2.0:0.1:2.0 1e20],'-','color',[0.6 0.6 0.6]);
% $$$ contour(X,Y,SSHr*100,[-1e20 -10.0:1:10.0 1e20],'-','color',[0.6 0.6 0.6]);
% $$$ plot(xreg,yreg,'-g','LineWidth',2);
% $$$ ylabel({'   Across-shelf ','Distance (km)'});
% $$$ cb = colorbar;
% $$$ %ylabel(cb,{'5-15 day MOM SSH','  anomaly (cm)'});
% $$$ caxis([-0.5 0.5]);
% $$$ set(gca,'xdir','reverse');
% $$$ xlim([X(1,x2) X(1,x1)]);
% $$$ set(gca,'color','k');
% $$$ %text(X(1,x1)+35,350,'(b)','BackgroundColor','w','FontSize',12,'Margin',0.1);
% $$$ set(gca,'Position',poss(end-1,:));
% $$$ set(gca,'xticklabel',[]);
% $$$ 
% $$$ subplot(6,1,3);
% $$$ contourf(X,Y,ROMSBTr*100,[-1e20 -1:0.025:1 1e20],'-k','linestyle','none');
% $$$ hold on;
% $$$ contour(X,Y,ROMSSSHr*100,[-1e20 -2.0:0.1:2.0 1e20],'-','color',[0.6 0.6 0.6]);
% $$$ plot(xreg,yreg,'-g','LineWidth',2);
% $$$ ylabel({'   Across-shelf ','Distance (km)'});
% $$$ cb = colorbar;
% $$$ %ylabel(cb,{'5-15 day MOM SSH','  anomaly (cm)'});
% $$$ caxis([-0.5 0.5]);
% $$$ set(gca,'xdir','reverse');
% $$$ xlim([X(1,x2) X(1,x1)]);
% $$$ set(gca,'color','k');
% $$$ %text(X(1,x1)+35,350,'(b)','BackgroundColor','w','FontSize',12,'Margin',0.1);
% $$$ set(gca,'Position',poss(end-2,:));
% $$$ set(gca,'xticklabel',[]);
% $$$ 
% $$$ % $$$ subplot(6,1,3);
% $$$ % $$$ contourf(X,Y,ROMSSSHr*100,[-1e20 -1.2:0.05:0 1e20],'-k','linestyle','none');
% $$$ % $$$ hold on;
% $$$ % $$$ plot(xreg,yreg,'-g','LineWidth',2);
% $$$ % $$$ ylabel({'   Across-shelf ','Distance (km)'});
% $$$ % $$$ cb = colorbar;
% $$$ % $$$ %ylabel(cb,{'5-15 day BT SSH','  anomaly (cm)'});
% $$$ % $$$ caxis([-0.01 0]*100);
% $$$ % $$$ set(gca,'xdir','reverse');
% $$$ % $$$ xlim([X(1,x2) X(1,x1)]);
% $$$ % $$$ set(gca,'color','k');
% $$$ % $$$ %text(X(1,x1)+35,350,'(c)','BackgroundColor','w','FontSize',12,'Margin',0.1);
% $$$ % $$$ set(gca,'Position',poss(end-2,:));
% $$$ % $$$ set(gca,'xticklabel',[]);
% $$$ % $$$ 
% $$$ % $$$ %ROMS BAROTROPIC VELOCITY:
% $$$ % $$$ subplot(6,1,4);
% $$$ % $$$ contourf(X,Y,ROMSBTr*100,[-1e20 -1:0.05:1 1e20],'-k','linestyle','none');
% $$$ % $$$ hold on;
% $$$ % $$$ plot(xreg,yreg,'-g','LineWidth',2);
% $$$ % $$$ ylabel({'   Across-shelf ','Distance (km)'});
% $$$ % $$$ cb = colorbar;
% $$$ % $$$ %ylabel(cb,{'5-15 day BT SSH','  anomaly (cm)'});
% $$$ % $$$ caxis([-0.5 0.5]);
% $$$ % $$$ set(gca,'xdir','reverse');
% $$$ % $$$ xlim([X(1,x2) X(1,x1)]);
% $$$ % $$$ set(gca,'color','k');
% $$$ % $$$ %text(X(1,x1)+35,350,'(c)','BackgroundColor','w','FontSize',12,'Margin',0.1);
% $$$ % $$$ set(gca,'Position',poss(end-2,:));
% $$$ % $$$ set(gca,'xticklabel',[]);
% $$$ 
% $$$ subplot(6,1,4);
% $$$ contourf(X,Y,Tr,[-1e20 -2:0.25:4 1e20],'linestyle','none');
% $$$ hold on;
% $$$ plot(xreg,yreg,'-g','LineWidth',2);
% $$$ ylabel({'   Across-shelf ','Distance (km)'});
% $$$ cb = colorbar;
% $$$ if (depthi == depthf)
% $$$     %ylabel(cb,['Control $\theta$ at ' num2str(depthi) 'm ($^\circ$C)']);
% $$$ else
% $$$     %ylabel(cb,['Control $\theta$ ' num2str(depthi) 'm - ' num2str(depthf) ...
% $$$     %               'm ($^\circ$C)']);
% $$$ end
% $$$ caxis([-2 4]);
% $$$ set(gca,'xdir','reverse');
% $$$ xlim([X(1,x2) X(1,x1)]);
% $$$ set(gca,'color','k');
% $$$ %text(X(1,x1)+35,350,'(d)','BackgroundColor','w','FontSize',12,'Margin',0.1);
% $$$ set(gca,'Position',poss(end-3,:));
% $$$ set(gca,'xticklabel',[]);
% $$$ 
% $$$ subplot(6,1,5);
% $$$ contourf(X,Y,TSOSEr,[-1e20 -2:0.25:4 1e20],'linestyle','none');
% $$$ hold on;
% $$$ plot(xreg,yreg,'-g','LineWidth',2);
% $$$ ylabel({'   Across-shelf ','Distance (km)'});
% $$$ cb = colorbar;
% $$$ if (depthi == depthf)
% $$$     %ylabel(cb,['SOSE $\theta$ at ' num2str(depthi) 'm ($^\circ$C)']);
% $$$ else
% $$$     %ylabel(cb,['SOSE $\theta$ ' num2str(depthi) 'm - ' num2str(depthf) ...
% $$$     %               'm ($^\circ$C)']);
% $$$ end
% $$$ caxis([-2 4]);
% $$$ set(gca,'xdir','reverse');
% $$$ xlim([X(1,x2) X(1,x1)]);
% $$$ set(gca,'color','k');
% $$$ %text(X(1,x1)+35,350,'(e)','BackgroundColor','w','FontSize',12,'Margin',0.1);
% $$$ set(gca,'Position',poss(end-4,:));
% $$$ set(gca,'xticklabel',[]);
% $$$ 
% $$$ subplot(6,1,6);
% $$$ contourf(X,Y,TDr,[-1e20 -1:0.1:1 1e20],'linestyle','none');
% $$$ hold on;
% $$$ plot(xreg,yreg,'-g','LineWidth',2);
% $$$ ylabel({'   Across-shelf ','Distance (km)'});
% $$$ cb = colorbar;
% $$$ if (depthi == depthf)
% $$$     %ylabel(cb,['$\theta$ anomaly at ' num2str(depthi) 'm ($^\circ$C)']);
% $$$ else
% $$$     %ylabel(cb,{['$\theta$ anomaly ' num2str(depthi) 'm - ' num2str(depthf) ...
% $$$     %               'm'],['days ' num2str(dayi) '-' num2str(dayf) '($^\circ$C)']});
% $$$ end
% $$$ caxis([-1 1]);
% $$$ set(gca,'xdir','reverse');
% $$$ xlim([X(1,x2) X(1,x1)]);
% $$$ set(gca,'color','k');
% $$$ %text(X(1,x1)+35,350,'(f)','BackgroundColor','w','FontSize',12,'Margin',0.1);
% $$$ set(gca,'Position',poss(end-5,:));
% $$$ xlabel('Distance along 1000m isobath (km)');
% $$$ 
% $$$ % $$$ subplot(4,1,4);
% $$$ % $$$ contourf(X,Y,N2r,[-1e20 0:0.25e-6:5e-6 1e20],'linestyle','none');
% $$$ % $$$ hold on;
% $$$ % $$$ plot(xreg,yreg,'-g','LineWidth',2);
% $$$ % $$$ ylabel({'     Distance from','1000m isobath (km)   '});
% $$$ % $$$ xlabel('Distance along 1000m isobath (km)');
% $$$ % $$$ cb = colorbar;
% $$$ % $$$ if (depthi == depthf)
% $$$ % $$$     ylabel(cb,{'Control $N^2$ at  ',['   ' num2str(depthi) 'm (s$^{-2}$)']});
% $$$ % $$$ else
% $$$ % $$$     ylabel(cb,{'Control $N^2$ at  ',['   ' num2str(depthi) 'm - ' num2str(depthf) ...
% $$$ % $$$                'm ($^\circ$C)']});
% $$$ % $$$ end
% $$$ % $$$ caxis([0 5e-6]);
% $$$ % $$$ xlim([X(1,x1) X(1,x2)]);
% $$$ % $$$ set(gca,'color','k');
% $$$ % $$$ text(X(1,x1)+35,350,'(d)','BackgroundColor','w','FontSize',12,'Margin',0.1);
% $$$ 
% $$$ colormap(redblue(40));
% $$$ 
% $$$ 
% $$$ 
% $$$ %%% Along Isobath SOSE data:
% $$$ theta_file = [base 'sose_theta.ltmm.nc'];
% $$$ salt_file = [base 'sose_salt.ltmm.nc'];
% $$$ sose_lon = ncread(theta_file,'longitude')-360;
% $$$ sose_lat = ncread(theta_file,'latitude');
% $$$ sose_depth = ncread(theta_file,'depth');
% $$$ sose_time = ncread(theta_file,'time');
% $$$ 
% $$$ [tmp yL] = min(abs(sose_lat+58));
% $$$ sose_lat = sose_lat(1:yL);
% $$$ [x_rho,y_rho] = meshgrid(sose_lon,sose_lat);
% $$$ xL = length(sose_lon);
% $$$ tL = length(sose_time);
% $$$ tL = 4; %FIRST 4 months ~~ first 125 days.
% $$$ 
% $$$ %[tmp zL] = min(abs(sose_depth+2500));
% $$$ %sose_depth = sose_depth(1:zL);
% $$$ 
% $$$ dep = 670;
% $$$ [tmp di] = min(abs(sose_depth+dep))
% $$$ 
% $$$ sose_theta = mean(ncread(theta_file,'temperature',[1 1 di-1 1],[xL yL 3 tL]),4);
% $$$ sose_salt = mean(ncread(salt_file,'salt',[1 1 di-1 1],[xL yL 3 tL]),4);
% $$$ 
% $$$ rho = sw_dens0(sose_salt,sose_theta);
% $$$ 
% $$$ N2 = avg(-9.81/1025*diff(rho,[],3)./ ...
% $$$          diff(repmat(permute(sose_depth(di-1:di+1),[3 2 1]),[xL yL 1]),[],3),3);
% $$$ 
% $$$ temp = sose_theta(:,:,2);
% $$$ 
% $$$ %Get lon-lat rotated:
% $$$ Tr = zeros(Nw+1,sL);
% $$$ N2r = zeros(Nw+1,sL);
% $$$ 
% $$$ W = 400;
% $$$ for k=1:sL
% $$$     k
% $$$     [lonr, latr, Corners,cc,lc] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,k));    
% $$$ 
% $$$     Trot = interp2(x_rho,y_rho,temp',lonr,latr,'linear');
% $$$     Tr(:,k) = nanmean(Trot,2);
% $$$     N2rot = interp2(x_rho,y_rho,N2',lonr,latr,'linear');
% $$$     N2r(:,k) = nanmean(N2rot,2);
% $$$ end
% $$$ 
% $$$ figure;
% $$$ set(gcf,'Position',[529         168        1035         820]);
% $$$ set(gcf,'defaulttextfontsize',12);
% $$$ set(gcf,'defaultaxesfontsize',12);
% $$$ set(gcf,'color','w');
% $$$ 
% $$$ %xpts:
% $$$ [tmp x1] = min(abs(cseg(1,:)+240));
% $$$ [tmp x2] = min(abs(cseg(1,:))-0);
% $$$ 
% $$$ [X,Y] = meshgrid((1:sL)*-dx,cc);
% $$$ X = X + 4000;
% $$$ 
% $$$ W = 210;
% $$$ %Regions:
% $$$ xreg = zeros(5,length(SECS));
% $$$ yreg = zeros(5,length(SECS));
% $$$ for i=1:length(SECS)
% $$$     xreg(:,i) = [X(1,SECS{i}(1)); ...
% $$$                  X(1,SECS{i}(end)); ...
% $$$                  X(1,SECS{i}(end)); ...
% $$$                  X(1,SECS{i}(1)); ...
% $$$                  X(1,SECS{i}(1))];
% $$$     yreg(:,i) = [-Wm; -Wm; W; W; -Wm];
% $$$ end
% $$$ 
% $$$ subplot(4,1,3);
% $$$ contourf(X,Y,Tr,[-1e20 -2:0.25:4 1e20],'linestyle','none');
% $$$ hold on;
% $$$ plot(xreg,yreg,'-g','LineWidth',2);
% $$$ ylabel({'     Distance from','1000m isobath (km)   '});
% $$$ cb = colorbar;
% $$$ ylabel(cb,'Control $\theta$ at 620m ($^\circ$C)');
% $$$ caxis([-2 3]);
% $$$ xlim([X(1,x1) X(1,x2)]);
% $$$ set(gca,'color','k');
% $$$ 
% $$$ subplot(4,1,4);
% $$$ contourf(X,Y,N2r,[-1e20 0:0.25e-6:3e-6 1e20],'linestyle','none');
% $$$ hold on;
% $$$ plot(xreg,yreg,'-g','LineWidth',2);
% $$$ ylabel({'     Distance from','1000m isobath (km)   '});
% $$$ xlabel('Distance along 1000m isobath (km)');
% $$$ cb = colorbar;
% $$$ ylabel(cb,{'Control $N^2$ at  ','   620m (s$^{-2}$)'});
% $$$ caxis([0 3e-6]);
% $$$ xlim([X(1,x1) X(1,x2)]);
% $$$ set(gca,'color','k');
% $$$ colormap(redblue(40));
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ %%%%%%%%%%%%%%OLD PLOTTING: 
% $$$ 
% $$$ %% 3 section plots movie:
% $$$ for dayi=5:5:100
% $$$     dayi = 150; %start day (actual)
% $$$     [tmp tii] = min(abs(time-dayi));
% $$$     dayf = 155; %end day (actual)
% $$$     [tmp tif] = min(abs(time-(dayf-5)));
% $$$     
% $$$     ti = [tii tif];
% $$$     
% $$$     %Get and plot variables along segments:
% $$$     files = {C_file, W_file, D_file};
% $$$     segs = {WANseg,OTseg,OT2seg};
% $$$     
% $$$     [lonr, latr, Corners,ccd,lcd] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,segs{1}(1)));    
% $$$     
% $$$     temp = zeros(length(ccd),zL,length(files),length(segs));
% $$$     dens = zeros(length(ccd),zL,length(files),length(segs));
% $$$     cc = zeros(length(ccd),zL,length(files),length(segs));
% $$$     lc = zeros(length(ccd),zL,length(files),length(segs));
% $$$     SSH  = zeros(length(ccd),1,length(files),length(segs));
% $$$     dep = zeros(length(ccd),zL,length(segs));
% $$$     
% $$$     for kk=1:length(segs)
% $$$         seg = segs{kk};
% $$$         [lonr, latr, Corners,ccd,lcd] = get_lonlat_rotated(W,Wm,Nw,Nl,L,cseg(:,seg(1)));    
% $$$ 
% $$$         cnt = 0;
% $$$         for k=1:length(seg)
% $$$             sprintf(['Doing segment %02d of ' names{kk}],k)
% $$$             for i = 1:length(files)
% $$$                 [T, tmp] = get_rotated_field(files{i},'temp','',ti,lonr,latr,'t','3D',cseg(3,seg(k)));
% $$$                 temp(:,:,i,kk) = temp(:,:,i,kk) + squeeze(nanmean(T,2));
% $$$ 
% $$$                 [T, tmp] = get_rotated_field(files{i},'pot_rho_0','',ti,lonr,latr,'t','3D',cseg(3,seg(k)));
% $$$                 dens(:,:,i,kk) = dens(:,:,i,kk) + squeeze(nanmean(T,2));
% $$$ 
% $$$ 
% $$$                 [U, V] = get_rotated_field(files{i},'u','v',ti,lonr,latr,'c','3D',cseg(3,k));
% $$$ % $$$                 [U, V] = get_rotated_field(files{i},'press_force_u','press_force_v',ti,lonr,latr,'c','3D',cseg(3,k));
% $$$                 cc(:,:,i,kk) = cc(:,:,i,kk) + squeeze(nanmean(U,2));
% $$$                 lc(:,:,i,kk) = lc(:,:,i,kk) + squeeze(nanmean(V,2));
% $$$ 
% $$$                 [ssh, tmp] = get_rotated_field(files{i},'sea_level','',ti,lonr,latr,'t','2D',cseg(3,k));
% $$$                 SSH(:,:,i,kk) = SSH(:,:,i,kk) + squeeze(nanmean(ssh,2));
% $$$             end
% $$$             [Z, tmp] = get_rotated_field(files{1},'dzt','',ti,lonr,latr,'t','3D',cseg(3,seg(k)));
% $$$             dep(:,:,kk) = dep(:,:,kk) + squeeze(nanmean(Z,2));
% $$$             cnt = cnt+1;
% $$$         end
% $$$         temp(:,:,:,kk) = temp(:,:,:,kk)/cnt;
% $$$         dens(:,:,:,kk) = dens(:,:,:,kk)/cnt;
% $$$         dep(:,:,kk) = dep(:,:,kk)/cnt;
% $$$         cc(:,:,:,kk) = cc(:,:,:,kk)/cnt;
% $$$         lc(:,:,:,kk) = lc(:,:,:,kk)/cnt;
% $$$         SSH(:,:,:,kk) = SSH(:,:,:,kk)/cnt;
% $$$     end
% $$$     
% $$$ % $$$     %Axes position (3x2)
% $$$ % $$$     posm = [0.0736    0.5587    0.25    0.3061; ...
% $$$ % $$$             0.3680    0.5587    0.25    0.3061; ...
% $$$ % $$$             0.6639    0.5587    0.25    0.3061; ...
% $$$ % $$$             0.0736    0.1153    0.25    0.3061; ...
% $$$ % $$$             0.3680    0.1153    0.25    0.3061; ...
% $$$ % $$$             0.6639    0.1153    0.25    0.3061];
% $$$ % $$$     poss = [0.0736    0.9    0.25    0.088; ...
% $$$ % $$$             0.3680    0.9    0.25    0.088; ...
% $$$ % $$$             0.6639    0.9    0.25    0.088; ...
% $$$ % $$$             0.0736    0.45    0.25    0.088; ...
% $$$ % $$$             0.3680    0.45    0.25    0.088; ...
% $$$ % $$$             0.6639    0.45    0.25    0.088];
% $$$ 
% $$$     %Axes positions (3x3):
% $$$     posm = [0.06    0.6978    0.25    0.2044; ...
% $$$             0.3680    0.6978    0.25    0.2044; ...
% $$$             0.68    0.6978    0.25    0.2044; ...
% $$$             0.06    0.375    0.25    0.2044; ...
% $$$             0.3680    0.375    0.25    0.2044; ...
% $$$             0.68    0.375    0.25    0.2044; ...
% $$$             0.06    0.045    0.25    0.2044; ...
% $$$             0.3680    0.045    0.25    0.2044; ...
% $$$             0.68    0.045    0.25    0.2044];
% $$$     poss = [0.06    0.9211    0.25    0.067; ...
% $$$             0.3680    0.9211    0.25    0.067; ...
% $$$             0.68    0.9211    0.25    0.067; ...
% $$$             0.06    0.61    0.25    0.067; ...
% $$$             0.3680    0.61    0.25    0.067; ...
% $$$             0.68    0.61   0.25    0.067; ...
% $$$             0.06    0.28    0.25    0.067; ...
% $$$             0.3680    0.28    0.25    0.067; ...
% $$$             0.68    0.28   0.25    0.067];
% $$$     txtx = -180;
% $$$     txty = -1250;
% $$$ 
% $$$     %figure;
% $$$     clf;
% $$$     %    set(gcf,'Position',[1          36        1920         970]);
% $$$     set(gcf,'Position',[217          40        1238         963]);
% $$$     set(gcf,'defaulttextfontsize',10);
% $$$     set(gcf,'defaultaxesfontsize',10);
% $$$     labels = {'CNTRL','WIND','DIFF'};
% $$$     xlims = [-200 210];
% $$$     xtic  = 100;
% $$$     ylims = [-1500 0];
% $$$     ytic = 500;
% $$$     colvar = temp;
% $$$     cntvar = dens;
% $$$     cntvar2 = lc;
% $$$     caxs = {[-2 4],[-2 4],[-0.5 0.5]}; %TEMP
% $$$     cint = {{[1020:0.04:1040],[400 400]},...
% $$$             {[1020:0.04:1040],[400 400]},...
% $$$             {[0.0025:0.0025:0.2],[-0.2:0.0025:-0.0025]}}; %isopycnals
% $$$     cint2 = {{[0.02:0.02:4],[-4:0.02:-0.02]},...
% $$$              {[0.02:0.02:4],[-4:0.02:-0.02]},...
% $$$              {[0.005:0.005:4],[-4:0.005:-0.005]}}; %along-coast velocity
% $$$     
% $$$ % $$$     colvar = cc;
% $$$ % $$$     cntvar = dens;
% $$$ % $$$     caxs = {[-3 3],[-3 3],[-1 1]}; %PR
% $$$ % $$$     cint = {{[1026.6:0.05:1027.8],[400 400]},...
% $$$ % $$$             {[1026.6:0.05:1027.8],[400 400]},...
% $$$ % $$$             {[0.005:0.005:0.06],[-0.06:0.005:-0.005]}}; %isopycnals
% $$$     
% $$$     SSHy_all = {{[-160 -130],[-160 -130],[-4.5 0.5]}, ...
% $$$                 {[-170 -140],[-170 -140],[-4.5 0.5]}, ...
% $$$                 {[-160 -130],[-160 -130],[-4.5 0.5]}};
% $$$     
% $$$     for kk=1:length(segs)
% $$$         
% $$$         SSHy = SSHy_all{kk};
% $$$         Z = -cumsum(dep(:,:,kk),2);
% $$$         X = repmat(ccd',[1 zL]);
% $$$         for i = 1:3
% $$$             ii = i;
% $$$             ii = i+3*(kk-1);
% $$$             ax = axes('Position',posm(ii,:));
% $$$             pcolPlot(X,Z,colvar(:,:,i,kk));
% $$$             hold on;
% $$$             contour(X,Z,cntvar(:,:,i,kk),cint{i}{1},'-k');
% $$$             contour(X,Z,cntvar(:,:,i,kk),cint{i}{2},'--k');
% $$$             contour(X,Z,cntvar2(:,:,i,kk),cint2{i}{1},'-m');
% $$$             contour(X,Z,cntvar2(:,:,i,kk),cint2{i}{2},'--m');
% $$$             caxis(caxs{i});
% $$$             ylim(ylims);
% $$$             xlim(xlims);
% $$$             
% $$$             text(txtx,txty,{names{kk}, ...
% $$$                           [labels{i} ' day ' num2str(dayi)]},'FontSize',10);
% $$$             set(gca,'ytick',[ylims(1):ytic:ylims(2)]);
% $$$             set(gca,'xtick',[xlims(1):xtic:xlims(2)]);
% $$$             set(gca,'xticklabel',[]);
% $$$             set(gca,'yticklabel',[]);
% $$$             if (ii == 1 | ii == 4 | ii == 7)
% $$$                 ylabel('Depth (m)');
% $$$                 set(gca,'yticklabel',[ylims(1):ytic:ylims(2)]);
% $$$             end
% $$$             if (ii >= 7)
% $$$                 set(gca,'xticklabel',[xlims(1):xtic:xlims(2)]);
% $$$                 xlabel('Distance from 1000m isobath (km)');    
% $$$             end    
% $$$             set(gca,'Position',posm(ii,:));
% $$$ 
% $$$             ax = axes('Position',poss(ii,:));
% $$$             plot(X(:,1),SSH(:,1,i,kk)*100,'-k');
% $$$             ylim(SSHy{i});
% $$$             xlim(xlims);
% $$$             set(gca,'xticklabel',[]);
% $$$             ylabel('SSH (cm)','FontSize',10);
% $$$         end
% $$$     end
% $$$ 
% $$$     name = sprintf('CC_movie/frame_dens%05d.png',ti);
% $$$     export_fig(name,'-painters');
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ 
