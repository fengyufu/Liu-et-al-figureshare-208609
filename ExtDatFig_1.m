clear
clc;

%% Loading data and converting the projection of data for figure plotting

load([pwd,'\input\data\target_reference_pixels'])
load([pwd,'\input\data\pixels_r_num']) % The number of reference pixels within a 1 degree window

NumMat = nan(3600,7200);
NumMat(MASKr==1) = 1;
NumMat(MASKt==1) = 2;

latlim = [40 90];
lonlim = [-180 180];

R = georasterref('RasterSize', [1000 7200], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B,RB] = geocrop(NumMat(1:1000,:),R,latlim,lonlim);
[lat,lon] = geographicGrid(RB);

%% Figure plotting 

cmap = [1,0,0; 0,1,0];

f = figure('Position',[0 0 800 600]);

% Figure a
ha = axes('Position',[0.13 0.11 0.49 0.82]);

load coastlines.mat
axesm('eqaazim','MapLatLimit',latlim,"MapLonLimit",lonlim)
axis off
framem on
gridm off
mlabel off
plabel off;

geoshow('landareas.shp','FaceColor',[0.96,0.96,0.96])
hold on

pcolorm(lat,lon,B)
colormap(gca,cmap)

% legend
axes('Position',[0.14 0.16 0.15 0.1]);
set(gca,'Color','none')
xlim([0,1])
ylim([0,1])
patch([0.3 0.5 0.5 0.3], [0.25 0.25 0.425 0.425], [1 0 0])
patch([0.3 0.5 0.5 0.3], [0.475 0.475 0.65 0.65], [0 1 0])
axis off

annotation(f,'textbox',[0.217 0.196 0.133 0.048],'String',{'Target pixels'},...
    'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.218 0.165 0.297 0.048],'String',{'Reference pixels (no FC change)'},...
    'LineStyle','none','FontSize',11);

% Figure b

mean = nanmean(rnum(:),"all"); % 35
ind35 = find(rnum == 35);
n = ind35(100);
[I,J] = ind2sub([3600,7200],find(MASKt==1));
MASKr_ck = MASKr(I(n)-9:I(n)+10,J(n)-10:J(n)+10);
MASKr_ck = double(MASKr_ck);
MASKr_ck(MASKr_ck==0) = nan;
MASKr_ck(10,10) = 2;

hb = axes('Position',[0.682 0.554 0.190 0.256]);
pcolor(MASKr_ck')
colormap(gca,cmap)
set(hb,'LineWidth',1,'YTick',[],'XTick',[])

% Figure c
hc = axes('Position',[0.682 0.214 0.192 0.286]);

histogram(rnum,'FaceAlpha',0.85,'EdgeColor','none','Normalization','probability')
xlim([0 150])
ylim([0 3e-4])
xlabel('Reference pixel Count')

set(hc,'FontSize',12,'box','off','LineWidth',1)