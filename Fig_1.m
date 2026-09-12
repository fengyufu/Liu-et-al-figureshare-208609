%% Fig. 1. Asymmetric sensitivity of ΔTAC to boreal forest gain and loss.

%% 
clear
clc;

%% Loading data

load([pwd,'\input\data\Tfcc_BF'])
load([pwd,'\input\data\Tfcc_BF_1d']) % ΔTAC in 1 degree for visualization
load([pwd,'\input\data\Koppen_BF'])
load([pwd,'\input\data\FC_diff_3y'])
FC_diff = FC_diff./100; % convert the range of forest cover fraction from 0-100(%) to 0-1 
load([pwd,'\input\data\alef_100']) % local sensitivity calculated in 1d x 1d window, 
                                   % TAC change per 100% forest cover change

%% Selecting pixels in boreal region

indnan = find(K~=1);

FC_diff(indnan) = nan;
TfccMap(indnan) = nan;
alef(indnan) = nan;

% classify local sensitivties into 12 categories
alef_cls = nan(3600,7200);
ind1 = find((K==1)&(alef<-0.4)&(p < 0.05)); alef_cls(ind1) = 1;
ind2 = find((K==1)&(alef>=-0.4)&(alef<-0.3)&(p < 0.05)); alef_cls(ind2) = 2;
ind3 = find((K==1)&(alef>=-0.3)&(alef<-0.2)&(p < 0.05)); alef_cls(ind3) = 3;
ind4 = find((K==1)&(alef>=-0.2)&(alef<-0.1)&(p < 0.05)); alef_cls(ind4) = 4;
ind5 = find((K==1)&(alef>=-0.1)&(alef<0)&(p < 0.05)); alef_cls(ind5) = 5;
ind6 = find((K==1)&(alef<0)&(p > 0.05)); alef_cls(ind6) = 6;
ind7 = find((K==1)&(alef>=0)&(p > 0.05)); alef_cls(ind7) = 7;
ind8 = find((K==1)&(alef>=0)&(alef<0.1)&(p < 0.05)); alef_cls(ind8) = 8;
ind9 = find((K==1)&(alef>=0.1)&(alef<0.2)&(p < 0.05)); alef_cls(ind9) = 9;
ind10 = find((K==1)&(alef>=0.2)&(alef<0.3)&(p < 0.05)); alef_cls(ind10) = 10;
ind11 = find((K==1)&(alef>=0.3)&(alef<0.4)&(p < 0.05)); alef_cls(ind11) = 11;
ind12 = find((K==1)&(alef>=0.4)&(p < 0.05)); alef_cls(ind12) = 12;

%% Converting the projection of data for figure plotting

latlim = [40 90];
lonlim = [-180 180];

R1 = georasterref('RasterSize', [1000 7200], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B1,RB1] = geocrop(TfccMap(1:1000,:),R1,latlim,lonlim);
[lat1,lon1] = geographicGrid(RB1);

R2 = georasterref('RasterSize', [1000 7200], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B2,RB2] = geocrop(FC_diff(1:1000,:),R2,latlim,lonlim);
[lat2,lon2] = geographicGrid(RB2);

R3 = georasterref('RasterSize', [1000 7200], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B3,RB3] = geocrop(alef_cls(1:1000,:),R3,latlim,lonlim);
[lat3,lon3] = geographicGrid(RB3);

% ensure the same region showing both for ΔTAC and ΔFC
indnan2 = find(~isfinite(B1));
B2(indnan2) = nan;

%% Colormap

load([pwd,'\input\cmaps\cmap1'])
load([pwd,'\input\cmaps\cmap2'])
load([pwd,'\input\cmaps\cmap3'])

%% Figure plotting

f = figure('Position',[0 0 1000 600]);

% Figure a
ha = axes('Position',[-0.0055 0.50 0.4 0.5]);

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
range = [-0.3,0.3];
clim(range)
cb1 = colorbar('southoutside');
set(cb1,'Position',[0.071 0.48 0.25 0.025],...
    'Ticks',[-0.3 -0.2 -0.1 0 0.1 0.2 0.3],...
    'TickLength',0.015,...
    'FontSize',11, ...
    'Label','ΔTAC');

% Figure b
hb = axes('Position',[0.31 0.50 0.4 0.5]);
hold(hb,'on');

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
colormap(gca,cmap2)
range = [-0.15,0.15];
clim(range)
cb2 = colorbar('southoutside');
set(cb2,'Position',[0.38 0.48 0.25 0.025],...
    'Ticks',[-0.15 -0.1 -0.05 0 0.05 0.1 0.15],...
    'TickLength',0.015,...
    'TickLabels',{'-15','-10','-5','0','5','10','15'},...
    'FontSize',11);

% Figure c
hc = axes('Position',[0.62 0.50 0.4 0.5]);
axis off
hold(hc,'on');

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
colormap(gca,cmap3)
cb3 = colorbar('southoutside');
set(cb3,'Position',[0.7 0.48 0.25 0.025],...
    'Ticks',[5.5835 7.4169],...
    'TickLength',0.055,...
    'TickLabels',{'0','0'},...
    'FontSize',11);

pro_neg = (length(find(alef<0))./length(find(isfinite(alef)))).*100;
pro_pos = (length(find(alef>0))./length(find(isfinite(alef)))).*100;

annotation(f,'textbox',[0.683 0.43 0.036 0.048],'String',{'-5'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.703 0.43 0.036 0.048],'String',{'-4'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.725 0.43 0.036 0.048],'String',{'-3'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.746 0.43 0.036 0.048],'String',{'-2'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.766 0.43 0.036 0.048],'String',{'-1'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.852 0.43 0.036 0.048],'String',{'1'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.873 0.43 0.036 0.048],'String',{'2'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.894 0.43 0.036 0.048],'String',{'3'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.914 0.43 0.036 0.048],'String',{'4'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.934 0.43 0.036 0.048],'String',{'5'},'LineStyle','none','FontSize',11);

annotation(f,'textbox',[0.84 0.5 0.05 0.05],'String',{'sig.'},'LineStyle','none','FontSize',12);
annotation(f,'textbox',[0.76 0.5 0.062 0.05],'String',{'insig.'},'LineStyle','none','FontSize',12);

annotation(f,'textbox',[0.94 0.47 0.058 0.09],'String','(×10^-^3)',...
    'LineStyle','none','FontSize',12,'FitBoxToText','off');

annotation(f,'textbox',[0.93 0.54 0.06 0.053],'Color',[0.85 0.33 0.10],...
    'String',[num2str(round(pro_pos)),'%'],'LineStyle','none','FontWeight','bold','FontSize',13);
annotation(f,'textbox',[0.93 0.58 0.06 0.05],'Color',[0.15 0.54 0.8],...
    'String',[num2str(round(pro_neg)),'%'],'LineStyle','none','FontWeight','bold','FontSize',13);

%% Figrue plotting bars

load([pwd,'\input\data\FC_multiyear_avg'])
FC_avg = FC_avg./100; % convert the range of forest cover fraction from 0-100(%) to 0-1 

IND = find((K==1)&(FC_avg > 0.05)); % inclusion of potential forested areas

FC_diff_b = nan(3600,7200);
FC_diff_b(IND) = FC_diff(IND);

edges = -0.25:0.02:0.25;
[N,edges,bin] = histcounts(FC_diff_b,edges);
X = -0.24:0.02:0.24;
Y = nan(length(edges)-1,1);
V = nan(length(edges)-1,1);
SEM = nan(length(edges)-1,1);
Y_ci = nan(length(edges)-1,2);
err = nan(length(edges)-1,2);
numsize = nan(length(edges)-1,1);

cellsize = 0.05;
lat(:,1) = -90+cellsize/2:cellsize:90-cellsize/2;
lon(1,:) = -180+cellsize/2:cellsize:180-cellsize/2;
[LON,LAT] = meshgrid(lon,flipud(lat));
W = cdtarea(LAT,LON,'km2');


for i = 1:length(edges)-1
    
    ind = find(bin==i);
    Y(i) = wmean(Tfcc(ind),W(ind),'omitnan');
    numsize(i) = length(find(isfinite(Tfcc(ind))));
    V(i) = (nanvar(Tfcc(ind),W(ind))^(0.5));
    SEM(i) = V(i)/sqrt(length(Tfcc(ind))); 
    ts = tinv([0.025  0.975],length(Tfcc(ind))-1); 
    err(i,:) = ts*SEM(i); % Y_ci(i,:) = Y(i) + ts*SEM(i);

end

% figrue d
hd = axes('Position',[0.10 0.095 0.55 0.28]);

X1 = X(1:13); X2 = X(13:end);
Y1 = Y(1:13); Y2 = Y(13:end);
b1 = bar(X1,Y1,0.8);

cmap4 = cat(1,cmap1(end-9,:),cmap1(10,:));

b1.FaceColor = cmap4(1,:);
b1.FaceAlpha = 0.5;
b1.EdgeColor = "none";
box on
hold on
b2 = bar(X2,Y2,0.8);
b2.FaceColor = cmap4(2,:);
b2.FaceAlpha = 0.5;
b2.EdgeColor = "none";
box on
hold on
er1 = errorbar(X1,Y1,err(1:13,1),err(1:13,2), ...
    'Linestyle','None','LineWidth',1,'CapSize',5,'Color',cmap4(1,:));
er2 = errorbar(X2,Y2,err(13:end,1),err(13:end,2), ...
    'Linestyle','None','LineWidth',1,'CapSize',5,'Color',cmap4(2,:));
ylim([-0.12, 0.12]);

plot([0,0],[-0.12,0.12],'k','Linewidth',0.7,'LineStyle','--','Marker','none')

set(hd,'FontSize',12,'LineWidth',1,'XTick',[-0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25], ...
    'XTickLabel',{'-25','-20','-15','-10','-5','0','5','10','15','20','25'}, ...
    'YTick', [-0.12 -0.06 0 0.06 0.12], ...
    'YTickLabel', {'-0.12','-0.06','0','0.06','0.12'});

ylabel({'\DeltaTAC'});
xlabel({'\DeltaFC (%)'});

% legend
annotation(f,'rectangle',...
    [0.54 0.33 0.025 0.02],'LineStyle','none','FaceColor',[0.70 0.59 0.47]);

annotation(f,'rectangle',...
    [0.54 0.30 0.025 0.02],'LineStyle','none','FaceColor',[0.50 0.68 0.65]);

annotation(f,'textbox',[0.56 0.285 0.10 0.05],'String',{'Forest gain'},'LineStyle','none',...
    'FontSize',12);

annotation(f,'textbox',[0.56 0.315 0.10 0.05],'String',{'Forest loss'},'LineStyle','none',...
    'FontSize',12);


% Figure e

indl = find((K==1)&(FC_avg>0.05)&(FC_diff_b<0));
indg = find((K==1)&(FC_avg>0.05)&(FC_diff_b>0));

dataX1 = FC_diff_b(indl);
dataY1 = Tfcc(indl);
[xData1, yData1] = prepareCurveData( dataX1, dataY1 );
[fitresult1, gof1, output1] = fit( xData1, yData1, 'poly1');
mdl1 = fitlm(xData1, yData1);
tValue1 = mdl1.Coefficients.tStat;

dataX2 = FC_diff_b(indg);
dataY2 = Tfcc(indg);
[xData2, yData2] = prepareCurveData( dataX2, dataY2 );
[fitresult2, gof2, output2] = fit( xData2, yData2, 'poly1');
mdl2 = fitlm(xData2, yData2);
tValue2 = mdl2.Coefficients.tStat;

cmap5 = [182,110,151; 124,134,65]./255;

hd = axes('Position',[0.75 0.10 0.22 0.28]);

x=[1,2]; 
y = [fitresult1.p1,fitresult2.p1];

err3 = [0.0139, 0.0125];

b3 = bar(x(1),y(1),0.5);
hold on 
er3 = errorbar(x(1),y(1),err3(1));
b3.FaceColor = cmap5(1,:);
b3.FaceAlpha = 0.7;
b3.EdgeColor = 'none';
b3.LineWidth = 1;
er3.LineWidth = 1.5;
er3.Color = cmap5(1,:);
er3.CapSize = 12;

b4 = bar(x(2),y(2),0.5);
hold on 
er4 = errorbar(x(2),y(2),err3(2));
b4.FaceColor = cmap5(2,:);
b4.FaceAlpha = 0.7;
b4.EdgeColor = 'none';
b4.LineWidth = 1.5;
er4.LineWidth = 1.5;
er4.Color = cmap5(2,:);
er4.CapSize = 12;

ylim([-0.6,0])
set(gca, 'YDir', 'reverse');

set(hd,'Color','none','FontSize',12,'LineWidth',1,'TickLength',...
    [0.02 0.025],'XTick',[1 2],'XTickLabel',{'',''},'YTick',[-0.6 -0.4 -0.2 0],...
    'YTickLabel',{'-6','-4','-2','0'},'Box','off'); 
annotation(f,'textbox',[0.75 0.337 0.074 0.06],'String',{'(×10^-^3)'},...
    'LineStyle','none','FontSize',12); % convertion into sensitivty of ΔTAC to 1% forest cover change (x10^-2)


ylabel({'Regional sensitivity ','of \DeltaTACto \DeltaFC'});

annotation(f,'textbox',[0.775 0.018 0.095 0.08],'String',{'Forest loss\newline    pixels'},...
    'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.855 0.017 0.095 0.08],'String',{'Forest gain\newline    pixels'},...
    'LineStyle','none','FontSize',11);

