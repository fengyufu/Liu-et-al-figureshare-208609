%% Supplementary Fig. 10. Spatial variability of climatic conditions and their variability within 1°x1° windows.

%% 
clear
clc;

%% Computing spatial variability of mean annual temperature

load([pwd,'\input\data\Koppen_BF_010d'])
load([pwd,'\input\data\t2m_ERA5Land'])


t2m_avg = nanmean(Xsum,3)- 273.15;
indnan = find(K_010==0);
t2m_avg(indnan) = nan;


rstart_pre = 1:10:1800;
rstart = floor(rstart_pre);

cstart_pre = 1:10:3600;
cstart = floor(cstart_pre);

var1 = nan(180,360);

% loop on moving window
for i = 1:length(rstart)
    for j = 1:length(cstart)

        chunk = t2m_avg(rstart(i):rstart(i)+9, cstart(j):cstart(j)+9);
        var1(i,j) = nanstd(chunk,0,'all');

    end
end

latlim = [40 90];
lonlim = [-180 180];

R1 = georasterref('RasterSize', [50 360], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B1,RB1] = geocrop(var1(1:50,:),R1,latlim,lonlim);
[lat1,lon1] = geographicGrid(RB1);


load([pwd,'\input\cmaps\cmap1'])

f = figure('position',[50,50,600,600]); 

% figure a
ha = axes('Position',[0.056 0.582 0.480 0.408]);

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
colormap(gca,cmap1)
clim([0,2])
clb = colorbar('southoutside');
set(clb,'Position',[0.122 0.582 0.35 0.02],'Ticks',[0 0.5 1 1.5 2],...
    'TickLabels',{'0','0.5','1.0','1.5','2.0'},'FontSize',10.8);

%% Computing spatial variability of temperature inter-annual variability

t2m_inter = nanstd(Xsum-273.15,0,3);
indnan = find(K_010==0);
t2m_inter(indnan) = nan;

rstart_pre = 1:10:1800;
rstart = floor(rstart_pre);

cstart_pre = 1:10:3600;
cstart = floor(cstart_pre);

var1 = nan(180,360);
% loop on moving window
for i = 1:length(rstart)
    for j = 1:length(cstart)

        chunk = t2m_inter(rstart(i):rstart(i)+9, cstart(j):cstart(j)+9);
        var2(i,j) = nanstd(chunk,0,'all');

    end
end

latlim = [40 90];
lonlim = [-180 180];

R2 = georasterref('RasterSize', [50 360], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B2,RB2] = geocrop(var2(1:50,:),R2,latlim,lonlim);
[lat2,lon2] = geographicGrid(RB2);


% figure b 
hb = axes('Position',[0.463 0.582 0.48 0.41]);

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
colormap(gca,cmap1)
clim([0,0.08])
clb = colorbar('southoutside');
set(clb,'Position',[0.527 0.582 0.35 0.02],'FontSize',10.8)

%% Computing spatial variability of mean annual precipitation

load([pwd,'\input\data\tp_ERA5Land'])

tp_avg = nanmean(Xsum,3);
indnan = find(K_010==0);
tp_avg(indnan) = nan;

rstart_pre = 1:10:1800;
rstart = floor(rstart_pre);

cstart_pre = 1:10:3600;
cstart = floor(cstart_pre);

var1 = nan(180,360);
% loop on moving window
for i = 1:length(rstart)
    for j = 1:length(cstart)

        chunk = tp_avg(rstart(i):rstart(i)+9, cstart(j):cstart(j)+9);
        var3(i,j) = nanstd(chunk,0,'all');

    end
end

latlim = [40 90];
lonlim = [-180 180];

R3 = georasterref('RasterSize', [50 360], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B3,RB3] = geocrop(var3(1:50,:),R3,latlim,lonlim);
[lat3,lon3] = geographicGrid(RB3);

% figure c
load([pwd,'\input\cmaps\cmap2'])

hc = axes('Position',[0.056 0.096 0.480 0.408]);

load coastlines.mat
axesm('eqaazim','MapLatLimit',latlim,"MapLonLimit",lonlim)
axis off
framem on
gridm off
mlabel off
plabel off;
geoshow('landareas.shp','FaceColor',[0.96,0.96,0.96])
hold on

pcolorm(lat3,lon3,B3)
colormap(gca,cmap2)
clim([0,350])
clb = colorbar('southoutside');
set(clb,'Position',[0.125 0.095 0.35 0.02],'Ticks',[0 70 140 210 280 350],'FontSize',10.8)

%% Computing spatial variability of precipitation inter-annual variability

tp_inter = nanstd(Xsum,0,3);
indnan = find(K_010==0);
tp_inter(indnan) = nan;

rstart_pre = 1:10:1800;
rstart = floor(rstart_pre);

cstart_pre = 1:10:3600;
cstart = floor(cstart_pre);

var1 = nan(180,360);
% loop on moving window
for i = 1:length(rstart)
    for j = 1:length(cstart)

        chunk = tp_inter(rstart(i):rstart(i)+9, cstart(j):cstart(j)+9);
        var4(i,j) = nanstd(chunk,0,'all');

    end
end

latlim = [40 90];
lonlim = [-180 180];

R4 = georasterref('RasterSize', [50 360], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B4,RB4] = geocrop(var4(1:50,:),R4,latlim,lonlim);
[lat4,lon4] = geographicGrid(RB4);

% figure d
hd = axes('Position',[0.463 0.096 0.480 0.408]);

load coastlines.mat
axesm('eqaazim','MapLatLimit',latlim,"MapLonLimit",lonlim)
axis off
framem on
gridm off
mlabel off
plabel off;
geoshow('landareas.shp','FaceColor',[0.96,0.96,0.96])
hold on

pcolorm(lat4,lon4,B4)
colormap(gca,cmap2)
clim([0,50])
clb = colorbar('southoutside');
set(clb,'Position',[0.529 0.095 0.350 0.02],'FontSize',10.8)
