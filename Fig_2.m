%% Fig. 2. Sensitivity of forest resilience to forest cover change across different forest types.

%% 
clear
clc;

%% Loading data

load([pwd,'\input\data\FML'])
load([pwd,'\input\data\Koppen_BF'])

indnan = find((K==0)|(FML == 53)); % Masking off areas in other biomes or with agroforestry 
FML(indnan) = nan;

%% Converting the projection of data for figure plotting

latlim = [40 90];
lonlim = [-180 180];

R1 = georasterref('RasterSize', [1000 7200], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B1,RB1] = geocrop(FML(1:1000,:),R1,latlim,lonlim);
[lat1,lon1] = geographicGrid(RB1);

%% Figure plotting 

f = figure('position',[50,50,1200,600]);

cmap1 = [94 166 156; 194 207 162; 164 121 158; 164 121 158]./255;

% Figure a
axes('Position',[0.005 0.5 0.31 0.47]);

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

% legend 
axes('position',[0.2 0.5 0.1,0.1])
axis off
set(gca,'Color','none')
xlim([0,1])
ylim([0,1])
patch([0.3 0.5 0.5 0.3], [0.25 0.25 0.425 0.425], [0.6431  0.4745 0.6196])
patch([0.3 0.5 0.5 0.3], [0.475 0.475 0.65 0.65], [0.7608  0.8118  0.6353])
patch([0.3 0.5 0.5 0.3], [0.7 0.7 0.875 0.875], [0.3686  0.6510  0.6118])

annotation(f,'textbox',[0.255 0.567 0.032 0.038],...
    'String',{'IF'},'LineStyle','none','FontSize',12);
annotation(f,'textbox',[0.252 0.529 0.038 0.051],...
    'String',{'RF'},'LineStyle','none','FontSize',12);
annotation(f,'textbox',[0.252 0.505 0.037 0.051],...
    'String',{'PF'},'LineStyle','none','FontSize',12);


% Figure b & c

load([pwd,'\input\data\long-term_TAC'])
load([pwd,'\input\data\alef_100']) % local sensitivity calculated in 1d x 1d window, TAC change per 100% forest cover change
load([pwd,'\input\data\FML_1d']) 
load([pwd,'\input\data\FC_diff_3y'])
FC_diff = FC_diff./100; % convert the range of forest cover fraction from 0-100(%) to 0-1 
load([pwd,'\input\data\Tfcc_BF'],'Tfcc')

% Creating grid
cellsize = 0.05;
lat(:,1) = -90+cellsize/2:cellsize:90-cellsize/2;
lon(1,:) = -180+cellsize/2:cellsize:180-cellsize/2;
[LON,LAT] = meshgrid(lon,flipud(lat));
W = cdtarea(LAT,LON,'km2');

ind_pri_005 = find(FML == 11);
ind_rege_005 = find(FML == 20);
ind_plan_005 = find((FML == 31)|(FML == 32));

data = nan(3,2); % row:forest types col:merics 
data(1,1) = wmean(ACF(ind_pri_005),W(ind_pri_005),'omitnan');
data(2,1) = wmean(ACF(ind_rege_005),W(ind_rege_005),'omitnan');
data(3,1) = wmean(ACF(ind_plan_005),W(ind_plan_005),'omitnan');

var1 = nan(3,2);
var1(1,1) = std(ACF(ind_pri_005),W(ind_pri_005),'omitnan');
var1(2,1) = std(ACF(ind_rege_005),W(ind_rege_005),'omitnan');
var1(3,1) = std(ACF(ind_plan_005),W(ind_plan_005),'omitnan');

dataX1 = FC_diff(ind_pri_005);
dataY1 = Tfcc(ind_pri_005);
[xData1, yData1] = prepareCurveData( dataX1, dataY1 );
[fitresult1, gof1] = fit( xData1, yData1, 'poly1' );
mdl1 = fitlm(xData1, yData1);
tValue1 = mdl1.Coefficients.tStat;

dataX2 = FC_diff(ind_rege_005);
dataY2 = Tfcc(ind_rege_005);
[xData2, yData2] = prepareCurveData( dataX2, dataY2 );
[fitresult2, gof2] = fit( xData2, yData2, 'poly1' );
mdl2 = fitlm(xData2, yData2);
tValue2 = mdl2.Coefficients.tStat;

dataX3 = FC_diff(ind_plan_005);
dataY3 = Tfcc(ind_plan_005);
[xData3, yData3] = prepareCurveData( dataX3, dataY3 );
[fitresult3, gof3] = fit( xData3, yData3, 'poly1' );
mdl3 = fitlm(xData3, yData3);
tValue3 = mdl3.Coefficients.tStat;

data(1,2) = fitresult1.p1;
data(2,2) = fitresult2.p1;
data(3,2) = fitresult3.p1;

var1(1,2) = 0.0092;
var1(2,2) = 0.0063;
var1(3,2) = 0.0191;

hb = axes('Position',[0.34 0.56 0.19 0.37]);

b1 = bar([1 2 3],data(:,1),'FaceAlpha',0.95,'FaceColor',[0.51 0.69 0.82],...
    'EdgeColor',[0.61 0.70 0.83],'BarWidth',0.6);
hold on
er1 = errorbar([1 2 3],data(:,1),var1(:,1),'LineStyle','none','LineWidth',1.5,...
    'Color',[0.25 0.30 0.53],'CapSize',8);
ylim([0, 0.8]);
set(hb,'FontSize',12,'LineWidth',1,'TickLength',[0.015 0.025],'XTick',...
    [1 2 3],'XTickLabel',{'IF','RF','PF'},'Box','off');
ylabel({'Long-term TAC'},'FontSize',13);

hc = axes('Position',[0.60 0.56 0.19 0.37]);
b2 = bar([1 2 3],data(:,2),'FaceColor',[0.78 0.68 0.77], ...
    'EdgeColor',[0.78 0.68 0.77],'BarWidth',0.6);
hold on
er2 = errorbar([1 2 3],data(:,2),var1(:,2),'LineStyle','none','LineWidth',1.5,...
    'Color',[0.44 0.41 0.56],'CapSize',9);
ylabel({'Regional sensitivity of','\DeltaTAC to \DeltaFC'});
ylim([-0.6, 0]); 
set(gca, 'YDir', 'reverse');
set(hc,'FontSize',12,'LineWidth',1,'TickLength',[0.015 0.025],'XTick',...
    [1 2 3],'XTickLabel',{'IF','RF','PF'},'YTick',[-0.6 -0.4 -0.2 0],...
    'YTickLabel',{'-6','-4','-2','0'},'Box','off');
annotation(f,'textbox',[0.60 0.9 0.06 0.06],...
    'String',{'(×10^-^3)'},'LineStyle','none','FontSize',12); % convertion into sensitivty of ΔTAC to 1% forest cover change (x10^-2)

% Figure d,e,f 

load([pwd,'\input\data\FC_diff_3y_1d']) % Forest cover change calculated at 1 degree scale

ind_pri_l = find((FML_t == 11)&(isfinite(alef))&(FC_diff_t < 0));
ind_pri_g = find((FML_t == 11)&(isfinite(alef))&(FC_diff_t > 0));

ind_rege_l = find((FML_t == 20)&(isfinite(alef))&(FC_diff_t < 0));
ind_rege_g = find((FML_t == 20)&(isfinite(alef))&(FC_diff_t > 0));

ind_plan_l = find(((FML_t == 31)|(FML_t == 32))&(isfinite(alef))&(FC_diff_t < 0));
ind_plan_g = find(((FML_t == 31)|(FML_t == 32))&(isfinite(alef))&(FC_diff_t > 0));

diff1 = wmean(alef(ind_pri_l),W(ind_pri_l),'omitnan') - wmean(alef(ind_pri_g),W(ind_pri_g),'omitnan');  [h1,p1] = ttest2(alef(ind_pri_l),alef(ind_pri_g),"Alpha",0.05);
diff2 = wmean(alef(ind_rege_l),W(ind_rege_l),'omitnan') - wmean(alef(ind_rege_g),W(ind_rege_g),'omitnan'); [h2,p2] = ttest2(alef(ind_rege_l),alef(ind_rege_g),"Alpha",0.05);
diff3 = wmean(alef(ind_plan_l),W(ind_plan_l),'omitnan') - wmean(alef(ind_plan_g),W(ind_plan_g),'omitnan'); [h3,p3] = ttest2(alef(ind_plan_l),alef(ind_plan_g),"Alpha",0.05);

% figure plotting
hd = axes('Position',[0.09 0.15 0.19 0.29]);
histogram(alef(ind_pri_l),'DisplayName','Forest loss','FaceAlpha',0.8,...
    'EdgeColor',[1 1 1],'FaceColor',[0.6 0.44 0.67],'Normalization','pdf','NumBins',80);
hold on
histogram(alef(ind_pri_g),'DisplayName','Forest gain','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.35 0.68 0.38],'Normalization','pdf','NumBins',250);
xlim([-1.5,1])
ylim([0 2])
ylabel({'Probability density'});
set(hd,'FontSize',12,'LineWidth',1,'XTickLabel',...
    {'-15','-10','-5','0','5','10'},'YTickLabel',{'0','0.5','1.0','1.5','2.0'},'Box','off');
title({'IF  -0.00065 (*)'});

he = axes('Position',[0.34 0.15 0.19 0.29]);
histogram(alef(ind_rege_l),'FaceAlpha',0.8,'EdgeColor',[1 1 1], ...
    'FaceColor',[0.6 0.44 0.67],'Normalization','pdf','NumBins',100);
hold on
histogram(alef(ind_rege_g),'FaceAlpha',0.55,'EdgeColor',[1 1 1],...
    'FaceColor',[0.35 0.68 0.38],'Normalization','pdf','NumBins',100);
xlim([-1.5,1])
ylim([0 2])
ylabel({'Probability density'});
xlabel({'Local sensitivity of \DeltaTAC to \DeltaFC'});
title({'RF   0.00003'});
set(he,'FontSize',12,'LineWidth',1,'XTickLabel',...
    {'-15','-10','-5','0','5','10'},'YTickLabel',{'0','0.5','1.0','1.5','2.0'},'Box','off');

hf = axes('Position',[0.6 0.15 0.19 0.29]);
histogram(alef(ind_plan_l),'DisplayName','Forest loss','FaceAlpha',0.8,'EdgeColor',[1 1 1],...
    'FaceColor',[0.6 0.44 0.67],'Normalization','pdf','NumBins',100);
hold on
histogram(alef(ind_plan_g),'DisplayName','Forest gain','FaceAlpha',0.55,'EdgeColor',[1 1 1],...
    'FaceColor',[0.35 0.68 0.38],'Normalization','pdf','NumBins',350);
xlim([-1.5,1])
ylim([0 2])
ylabel({'Probability density'});
title({'PF  0.00119 (*)'});
set(hf,'FontSize',12,'LineWidth',1,'XTickLabel',...
    {'-15','-10','-5','0','5','10'},'YTickLabel',{'0','0.5','1.0','1.5','2.0'},'Box','off');
annotation(f,'textbox',[0.79 0.14 0.06 0.06],...
    'String',{'(×10^-^3)'},'LineStyle','none','FontSize',12); % convertion into sensitivty of ΔTAC to 1% forest cover change (x10^-2)


