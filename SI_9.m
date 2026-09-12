%% Supplementary Fig. 9. Comparison of TAC trends using moving windows of different lengths.

%% 
clear
clc;

%% Loading data
load([pwd,'\input\data\t2m_ERA5Land_multiyear_avg'])
t2m = Xavg-273.15;

load([pwd,'\input\data\tp_ERA5Land_multiyear_avg'])
tp = Xavg;

load([pwd,'\input\data\FC_multiyear_avg'])
load([pwd,'\input\data\Koppen_BF'])
load([pwd,'\input\data\kNDVI_missing_proportion'])

clear Xavg

%% Setting

indnan = find((FC_avg < 5)|(K ~= 1)|(MDN > 0.5));
t2m(indnan) = nan;
tp(indnan) = nan;

cellsize = 0.05;
lat1(:,1) = -90+cellsize/2:cellsize:90-cellsize/2;
lon1(1,:) = -180+cellsize/2:cellsize:180-cellsize/2;
[LON,LAT] = meshgrid(lon1,flipud(lat1));
W = cdtarea(LAT,LON,'km2');

xlimset = [-10 10; 200 1200];

minX = xlimset(1,1);
maxX = xlimset(1,2);
minY = xlimset(2,1);
maxY = xlimset(2,2);

%% TAC trends using 3-year moving window

f = figure('Position',[0 0 746 617]);

% figure a
load([pwd,'\input\data\TAC_rl3_trend_TS'])
load([pwd,'\input\cmap\cmap'])

Xtrend(indnan) = nan;

latlim = [40 90];
lonlim = [-180 180];

R1 = georasterref('RasterSize', [1000 7200], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B1,RB1] = geocrop(Xtrend(1:1000,:),R1,latlim,lonlim);
[lat1,lon1] = geographicGrid(RB1);

ha = axes('Position',[0 0.52 0.5 0.38]);

load coastlines.mat
axesm('eqaazim','MapLatLimit',latlim,"MapLonLimit",lonlim)
axis off
framem on
gridm off
mlabel off
plabel off;

geoshow('landareas.shp','FaceColor',[0.96,0.96,0.96])
hold on

pcolorm(lat1,lon1,B1)
colormap(gca,cmap)
clim([-0.02,0.02])


xx = t2m;
yy = tp;
zz = Xtrend;
ww = W;

xx(xx<minX) = nan;
xx(xx>maxX) = nan;
yy(yy<minY) = nan;
yy(yy>maxY) = nan;

[INDnan] = find(isnan(xx)|isnan(yy));
xx(INDnan) = [];
yy(INDnan) = [];
zz(INDnan) = [];

targetSize = [40 40];

xxBin = round( (xx-minX)/(maxX-minX)*(targetSize(1)-1) ) +1;
yyBin = round( (yy-minY)/(maxY-minY)*(targetSize(2)-1) ) +1;
    

Znum = accumarray([xxBin(:),yyBin(:)],~isnan(zz),targetSize,@nansum,nan);

% Re-gridding
Xvet = minX:((maxX-minX)/(targetSize(1)-1)):maxX;
Yvet = minY:((maxY-minY)/(targetSize(2)-1)):maxY;
[Xmat, Ymat] = meshgrid(Xvet, Yvet);

HB = nan(40,40);
PB = nan(40,40);
Z = nan(40,40);

for ibin = 1:targetSize(1)
    for jbin=1:targetSize(2)
        
        INDbin = find((xxBin==ibin)&(yyBin==jbin));
        zt = zz(INDbin);
        wt = ww(INDbin);

        if (~isempty(zt))&&(length(zt)<2)
            
            Z(ibin,jbin) = zt;
            
        elseif length(zt)>1 

            Z(ibin,jbin) = wmean(zt,wt,'omitnan');
            
            [hB,pB] = ttest(zt, 0, 0.05,'both');
            
            HB(ibin,jbin) = hB;
            PB(ibin,jbin) = pB;
        end
    end
end

hb = axes('Position',[0 0.1 0.5 0.38]);

Z(Znum<50) = nan; 
pcolor(Xmat,Ymat,Z'); shading flat; axis square;
axis([min(Xmat(:)) max(Xmat(:)) min(Ymat(:)) max(Ymat(:))])
set(gca,'layer','top')
colormap(gca,cmap)
clim([-0.02, 0.02])
xlabel('Mean annual temperature (℃)','Fontsize',12);
ylabel('Mean annual precipitation (mm)','Fontsize',12);
    
% Overlay significance
IND = find((HB'==1)&(isfinite(Z')));
Xind = Xmat(IND)+0.5*abs((Xmat(1,1)-Xmat(1,2)));
Yind = Ymat(IND)+0.5*abs((Ymat(1,1)-Ymat(2,1)));
hold on; plot(Xind, Yind,'.k','MarkerSize',4)

set(hb,'FontSize',12,'Layer','top','LineWidth',1)

clb = colorbar('southoutside');
set(clb,'Position',[0.372 0.524 0.265 0.017],'TickLength',0.02,'FontSize',10.8);

%% TAC trends using 5-year moving window

% figure b
load([pwd,'\input\data\TAC_rl5_trend_TS'])

Xtrend(indnan) = nan;

latlim = [40 90];
lonlim = [-180 180];

R2 = georasterref('RasterSize', [1000 7200], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B2,RB2] = geocrop(Xtrend(1:1000,:),R2,latlim,lonlim);
[lat2,lon2] = geographicGrid(RB2);

hc = axes('Position',[0.5 0.52 0.5 0.38]);

load coastlines.mat
axesm('eqaazim','MapLatLimit',latlim,"MapLonLimit",lonlim)
axis off
framem on
gridm off
mlabel off
plabel off;

geoshow('landareas.shp','FaceColor',[0.96,0.96,0.96])
hold on

pcolorm(lat2,lon2,B2)
colormap(gca,cmap)
clim([-0.02,0.02])


xx = t2m;
yy = tp;
zz = Xtrend;
ww = W;

xx(xx<minX) = nan;
xx(xx>maxX) = nan;
yy(yy<minY) = nan;
yy(yy>maxY) = nan;

[INDnan] = find(isnan(xx)|isnan(yy));
xx(INDnan) = [];
yy(INDnan) = [];
zz(INDnan) = [];

targetSize = [40 40];

xxBin = round( (xx-minX)/(maxX-minX)*(targetSize(1)-1) ) +1;
yyBin = round( (yy-minY)/(maxY-minY)*(targetSize(2)-1) ) +1;
    

Znum = accumarray([xxBin(:),yyBin(:)],~isnan(zz),targetSize,@nansum,nan);

% Re-gridding
Xvet = minX:((maxX-minX)/(targetSize(1)-1)):maxX;
Yvet = minY:((maxY-minY)/(targetSize(2)-1)):maxY;
[Xmat, Ymat] = meshgrid(Xvet, Yvet);

HB = nan(40,40);
PB = nan(40,40);
Z = nan(40,40);

for ibin = 1:targetSize(1)
    for jbin=1:targetSize(2)
        
        INDbin = find((xxBin==ibin)&(yyBin==jbin));
        zt = zz(INDbin);
        wt = ww(INDbin);

        if (~isempty(zt))&&(length(zt)<2)
            
            Z(ibin,jbin) = zt;
            
        elseif length(zt)>1 

            Z(ibin,jbin) = wmean(zt,wt,'omitnan');
            
            [hB,pB] = ttest(zt, 0, 0.05,'both');
            
            HB(ibin,jbin) = hB;
            PB(ibin,jbin) = pB;
        end
    end
end

hd = axes('Position',[0.5 0.1 0.5 0.38]);

Z(Znum<50) = nan; 
pcolor(Xmat,Ymat,Z'); shading flat; axis square;
axis([min(Xmat(:)) max(Xmat(:)) min(Ymat(:)) max(Ymat(:))])
set(gca,'layer','top')
colormap(gca,cmap)
clim([-0.02, 0.02])
xlabel('Mean annual temperature (℃)','Fontsize',12);
ylabel('Mean annual precipitation (mm)','Fontsize',12);
    
% Overlay significance
IND = find((HB'==1)&(isfinite(Z')));
Xind = Xmat(IND)+0.5*abs((Xmat(1,1)-Xmat(1,2)));
Yind = Ymat(IND)+0.5*abs((Ymat(1,1)-Ymat(2,1)));
hold on; plot(Xind, Yind,'.k','MarkerSize',4)

set(hd,'FontSize',12,'Layer','top','LineWidth',1)

