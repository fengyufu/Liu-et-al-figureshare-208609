%% Supplementary Fig. 5. Sensitivity of ΔTAC to ΔFC within pixels with distinct trajectories of forest gain.

%% Supplementary Fig. 5. Sensitivity of ΔTAC to ΔFC within pixels with distinct trajectories of forest gain. 

%% 

clear
clc;

%% Loading data 

load([pwd,'\input\data\FC_diff_3y_1d']) % Forest cover change calculated at 1 degree scale
load([pwd,'\input\data\Koppen_BF'])
load([pwd,'\input\data\alef_100']) % local sensitivity calculated in 1d x 1d window, 
                                   % TAC change per 100% forest cover change

indnan = find(K~=1); % Masking off areas in other biomes

FC_diff_t(indnan) = nan;
dif1_t(indnan) = nan; % Forest cover change in 2000-2010 period
dif2_t(indnan) = nan; % Forest cover change in 2011-2020 period

%% Classifying data

indg_gg = find((FC_diff_t > 0)&(dif1_t > 0)&(dif2_t > 0));
indg_lg = find((FC_diff_t > 0)&(dif1_t < 0)&(dif2_t > 0));
indg_gl = find((FC_diff_t > 0)&(dif1_t > 0)&(dif2_t < 0));

NumMat = nan(3600,7200);

NumMat(indg_gg) = 1;
NumMat(indg_lg) = 2;
NumMat(indg_gl) = 3;

%% Converting the projection of data for figure plotting

cellsize = 0.05;
lat(:,1) = -90+cellsize/2:cellsize:90-cellsize/2;
lon(1,:) = -180+cellsize/2:cellsize:180-cellsize/2;
[LON,LAT] = meshgrid(lon,flipud(lat));
W = cdtarea(LAT,LON,'km2');

latlim = [40 90];
lonlim = [-180 180];

R1 = georasterref('RasterSize', [1000 7200], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B1,RB1] = geocrop(NumMat(1:1000,:),R1,latlim,lonlim);
[lat1,lon1] = geographicGrid(RB1);

%% Figure plotting

f = figure('Position',[50,50,620,600]);

% Figure a
ha = axes('Position',[-0.033 0.506 0.505 0.487]);

cmap = [0.67 0.58 0.96; 0.14 0.82 0.88; 1 0.78 0];

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

% Figure b

data = nan(3,1);
data(1,1) = wmean(alef(indg_gg),W(indg_gg),'omitnan');
data(2,1) = wmean(alef(indg_lg),W(indg_lg),'omitnan');
data(3,1) = wmean(alef(indg_gl),W(indg_gl),'omitnan');

var = nan(3,1);
var(1,1) = std(alef(indg_gg),W(indg_gg),'omitnan');
var(2,1) = std(alef(indg_lg),W(indg_lg),'omitnan');
var(3,1) = std(alef(indg_gl),W(indg_gl),'omitnan');

SEM = nan(3,1);
SEM(1,1) = var(1,1)/sqrt(length(alef(indg_gg))); 
SEM(2,1) = var(2,1)/sqrt(length(alef(indg_lg))); 
SEM(3,1) = var(3,1)/sqrt(length(alef(indg_gl))); 

ts = nan(3,2);
ts(1,:) = tinv([0.025  0.975],length(alef(indg_gg))-1); 
ts(2,:) = tinv([0.025  0.975],length(alef(indg_lg))-1); 
ts(3,:) = tinv([0.025  0.975],length(alef(indg_gl))-1); 

err = ts.*SEM;

hb = axes('Position',[0.571 0.542 0.409 0.410]);

x=[1,2,3];

b1 = bar(x(1),data(1),'FaceAlpha',0.6,'EdgeColor',[1 1 1],'BarWidth',0.55,'FaceColor',[0.43,0.36,0.75]);
hold on 
er1 = errorbar(x(1),data(1),err(1,2),'LineWidth',1.8,'CapSize',8,'Color',[0.43,0.36,0.75]);

b2 = bar(x(2),data(2),'FaceAlpha',0.65,'EdgeColor',[1 1 1],'BarWidth',0.55,'FaceColor',[0.31,0.81,0.79]);
hold on 
er2 = errorbar(x(2),data(2),err(2,2),'LineWidth',1.8,'CapSize',8,'Color',[0.31,0.81,0.79]);

b3 = bar(x(3),data(3),'FaceAlpha',0.95,'EdgeColor',[1 1 1],'BarWidth',0.55,'FaceColor',[0.96,0.90,0.64]);
hold on 
er3 = errorbar(x(3),data(3),err(3,2),'LineWidth',1.8,'CapSize',8,'Color',[0.96,0.75,0.64]);

set(gca, 'YDir', 'reverse');
xlim([0,4])
ylim([-0.25,0])

ylabel({'Local sensitivity of \DeltaTAC to \DeltaFC'});
xlabel({''});

set(hb,'FontSize',12,'LineWidth',1,'TickLength',[0.005 0.025],'XTick',zeros(1,0), ...
    'Box','off','YTickLabel',{'-2.5','-2.0','-1.5','-1.0','-0.5','0'});

annotation(f,'textbox',[0.557 0.948 0.111 0.054],...
    'String',{'(×10^-^3)'},'LineStyle','none','FontSize',11); % convertion into sensitivty of ΔTAC to 1% forest cover change (x10^-2)