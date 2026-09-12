clc;
clear

%% ESA CCI 
% Loading data
load('D:\assignment\ForestResilience\response\ESA_CCI_LC\output\TfccTres')
load('D:\assignment\ForestResilience\response\ESA_CCI_LC\output\pfts_for_ts') % FC ranging from 0 to 1
FC_diff = nanmean(Xts(:,:,19:21),3) - nanmean(Xts(:,:,1:3),3);
load('D:\assignment\ForestResilience\koppen\output\Koppen_005')

IND = find((K==1)&(Xavg>0.05)); % Extract boreal forest 

FC_diff_b = nan(3600,7200);
FC_diff_b(IND) = FC_diff(IND);

edges = -0.25:0.02:0.25;
[N,edges,bin] = histcounts(FC_diff_b,edges);
X = -0.24:0.02:0.24;
Y = nan(length(edges)-1,1);
V = nan(length(edges)-1,1);
SEM = nan(length(edges)-1,1);
err = nan(length(edges)-1,2);

cellsize = 0.05;
lat(:,1) = -90+cellsize/2:cellsize:90-cellsize/2;
lon(1,:) = -180+cellsize/2:cellsize:180-cellsize/2;
[LON,LAT] = meshgrid(lon,flipud(lat));
W = cdtarea(LAT,LON,'km2');

for i = 1:length(edges)-1
    
    ind = find(bin==i);
    Y(i) = wmean(Tfcc(ind),W(ind),'omitnan');
    V(i) = (nanvar(Tfcc(ind),W(ind))^(0.5));
    SEM(i) = V(i)/sqrt(length(Tfcc(ind))); 
    ts = tinv([0.025  0.975],length(Tfcc(ind))-1); 
    err(i,:) = ts*SEM(i); 

end

indl = find((K==1)&(Xavg>0.05)&(FC_diff_b<0));
indg = find((K==1)&(Xavg>0.05)&(FC_diff_b>0));

dataX1 = FC_diff_b(indl);
dataY1 = Tfcc(indl);
[xData1, yData1] = prepareCurveData( dataX1, dataY1 );
[fitresult1, gof1] = fit( xData1, yData1, 'poly1' );

dataX2 = FC_diff_b(indg);
dataY2 = Tfcc(indg);
[xData2, yData2] = prepareCurveData( dataX2, dataY2 );
[fitresult2, gof2] = fit( xData2, yData2, 'poly1' );

% Figrue plotting 
f = figure('Position',[0 0 1300 832]);

% figure c,d
hc = axes('Position',[0.125 0.566 0.5 0.17]);

mycolor = [0.41	0.23 0.03; 0 0.29 0.25];
mycolor2 = [0.71 0.43 0.59; 0.49 0.53 0.25];

X1 = X(1:13); X2 = X(13:end);
Y1 = Y(1:13); Y2 = Y(13:end);
b1 = bar(X1,Y1,0.8);
b1.FaceColor = mycolor(1,:);
b1.FaceAlpha = 0.5;
b1.EdgeColor = "none";
box on
hold on
b2 = bar(X2,Y2,0.8);
b2.FaceColor = mycolor(2,:);
b2.FaceAlpha = 0.5;
b2.EdgeColor = "none";
box on
hold on
er1 = errorbar(X1,Y1,err(1:13,1),err(1:13,2), ...
    'Linestyle','None','LineWidth',1,'CapSize',5,'Color',mycolor(1,:));
er2 = errorbar(X2,Y2,err(13:end,1),err(13:end,2), ...
    'Linestyle','None','LineWidth',1,'CapSize',5,'Color',mycolor(2,:));
ylim([-0.2, 0.2]);

plot([0,0],[-0.2,0.2],'k','Linewidth',0.5,'LineStyle','--')

set(hc,'FontSize',14,'GridLineWidth',0.7,'LineWidth',1,'MinorGridAlpha',...
    0.2,'XMinorGrid','on','XTick',...
    [-0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25],'XTickLabel',...
    {'-25','-20','-15','-10','-5','0','5','10','15','20','25'},'YMinorGrid',...
    'on','YTick',[-0.2 -0.1 0 0.1 0.2],'YTickLabel',...
    {'-0.2','-0.1','0','0.1','0.2'});

ylabel({'\DeltaTAC of ','kNDVI'});
xlabel({'\DeltaFC derived from ESA CCI (%)'},'FontSize',14);


hd = axes('Position',[0.707,0.566,0.16,0.17]);

x=[1,2];   
y=[fitresult1.p1,fitresult2.p1];

err1 = [0.024, 0.024];
b3 = bar(x(1),y(1),0.5);
hold on 
er3 = errorbar(x(1),y(1),err1(1));
b3.FaceColor = mycolor2(1,:);
b3.FaceAlpha = 0.7;
b3.EdgeColor = 'none';
b3.LineWidth = 1;
er3.LineWidth = 1.5;
er3.Color = mycolor2(1,:);
er3.CapSize = 10;

b4 = bar(x(2),y(2),0.5);
hold on 
er4 = errorbar(x(2),y(2),err1(2));
b4.FaceColor = mycolor2(2,:);
b4.FaceAlpha = 0.7;
b4.EdgeColor = 'none';
b4.LineWidth = 1;
er4.LineWidth = 1.5;
er4.Color = mycolor2(2,:);
er4.CapSize = 10;

ylim([-0.8, 0]);
set(gca, 'YDir', 'reverse','Box','off','FontSize',13,'LineWidth',1,'XTick',[1 2],...
    'XTickLabel',{'',''},'YTick',[-0.8 -0.6 -0.4 -0.2 0],'YTickLabel',...
    {'-8','-6','-4','-2','0'});
ylabel({'Regional sensitivity ','of \DeltaTAC to \DeltaFC'},'FontSize',13);

annotation(f,'textbox',[0.706 0.707 0.05 0.04],...
    'String','(×10^-^3)','LineStyle','none','FontSize',12);

annotation(f,'textbox',[0.72 0.51 0.07 0.06],...
    'String',{'Forest loss\newline    pixels'},'LineStyle','none','FontSize',11);

annotation(f,'textbox',[0.78 0.51 0.07 0.06],...
    'String',{'Forest gain\newline    pixels'},'LineStyle','none','FontSize',11);

%% Variance
% Loading data
load('D:\assignment\ForestResilience\response\variance\output\CV_rl3_fcc.mat')
load('D:\assignment\ForestResilience\vcf\output\whole\FC_avg.mat') % FC ranging from 0 to 100
load('D:\assignment\ForestResilience\vcf\output\whole\FC_diff_3y.mat')
FC_diff = FC_diff./100;
load('D:\assignment\ForestResilience\koppen\output\Koppen_005.mat')

IND = find((K==1)&(FC_avg>5));

FC_diff_b = nan(3600,7200);
FC_diff_b(IND) = FC_diff(IND);

edges = -0.25:0.02:0.25;
[N,edges,bin] = histcounts(FC_diff_b,edges);
X = -0.24:0.02:0.24;
Y = nan(length(edges)-1,1);
V = nan(length(edges)-1,1);
SEM = nan(length(edges)-1,1);
err = nan(length(edges)-1,2);

for i = 1:length(edges)-1
    
    ind = find(bin==i);
    Y(i) = wmean(Cfcc(ind),W(ind),'omitnan');
    V(i) = (nanvar(Cfcc(ind),W(ind))^(0.5));
    SEM(i) = V(i)/sqrt(length(Cfcc(ind))); 
    ts = tinv([0.025  0.975],length(Cfcc(ind))-1); 
    err(i,:) = ts*SEM(i); % Y_ci(i,:) = Y(i) + ts*SEM(i);

end

indl = find((K==1)&(FC_avg>0.05)&(FC_diff_b<0));
indg = find((K==1)&(FC_avg>0.05)&(FC_diff_b>0));

dataX1 = FC_diff_b(indl);
dataY1 = Cfcc(indl);
[xData1, yData1] = prepareCurveData( dataX1, dataY1 );
[fitresult1, gof1] = fit( xData1, yData1, 'poly1' );

dataX2 = FC_diff_b(indg);
dataY2 = Cfcc(indg);
[xData2, yData2] = prepareCurveData( dataX2, dataY2 );
[fitresult2, gof2] = fit( xData2, yData2, 'poly1' );

% figure e,f
he = axes('Position',[0.125,0.319,0.50,0.17]);

X1 = X(1:13); X2 = X(13:end);
Y1 = Y(1:13); Y2 = Y(13:end);
b1 = bar(X1,Y1,0.8);
b1.FaceColor = mycolor(1,:);
b1.FaceAlpha = 0.5;
b1.EdgeColor = "none";
box on
hold on
b2 = bar(X2,Y2,0.8);
b2.FaceColor = mycolor(2,:);
b2.FaceAlpha = 0.5;
b2.EdgeColor = "none";
box on
hold on
er1 = errorbar(X1,Y1,err(1:13,1),err(1:13,2), ...
    'Linestyle','None','LineWidth',1,'CapSize',5,'Color',mycolor(1,:));
er2 = errorbar(X2,Y2,err(13:end,1),err(13:end,2), ...
    'Linestyle','None','LineWidth',1,'CapSize',5,'Color',mycolor(2,:));
ylim([-0.015, 0.015]);

plot([0,0],[-0.015,0.015],'k','Linewidth',0.5,'LineStyle','--')

set(he,'FontSize',14,'GridLineWidth',0.7,'LineWidth',1,'MinorGridAlpha',...
    0.2,'XMinorGrid','on','XTick',...
    [-0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25],'XTickLabel',...
    {'-25','-20','-15','-10','-5','0','5','10','15','20','25'},'YMinorGrid',...
    'on','YTick',[-0.015 -0.0075 0 0.0075 0.015],'YTickLabel',...
    {'-15.0','-7.5','0','7.5','15.0'});

ylabel({'\DeltaCV of ','kNDVI'});
xlabel({'\DeltaFC derived from MOD44B (%)'},'FontSize',14);

annotation(f,'textbox',[0.12 0.485 0.06 0.04],'String',{'(x10^-^3)'},...
    'LineStyle','none','FontSize',12);

hf = axes('Position',[0.7074,0.3188,0.16,0.17]);

x=[1,2];   
y=[fitresult1.p1,fitresult2.p1];

err1 = [0.0042, 0.071];
b3 = bar(x(1),y(1),0.5);
hold on 
er3 = errorbar(x(1),y(1),err1(1));
b3.FaceColor = mycolor2(1,:);
b3.FaceAlpha = 0.7;
b3.EdgeColor = 'none';
b3.LineWidth = 1;
er3.LineWidth = 1.5;
er3.Color = mycolor2(1,:);
er3.CapSize = 10;

b4 = bar(x(2),y(2),0.5);
hold on 
er4 = errorbar(x(2),y(2),err1(2));
b4.FaceColor = mycolor2(2,:);
b4.FaceAlpha = 0.7;
b4.EdgeColor = 'none';
b4.LineWidth = 1;
er4.LineWidth = 1.5;
er4.Color = mycolor2(2,:);
er4.CapSize = 10;

ylim([-0.08, 0]);
set(gca, 'YDir', 'reverse','Box','off','FontSize',13,'LineWidth',1,'XTick',[1 2],...
    'XTickLabel',{'',''},'YTick',[-0.08 -0.06 -0.04 -0.02 0],'YTickLabel',...
    {'-8','-6','-4','-2','0'});
ylabel({'Regional sensitivity ','of \DeltaTAC to \DeltaFC'},'FontSize',13);

annotation(f,'textbox',[0.706 0.459 0.05 0.04],...
    'String','(×10^-^4)','LineStyle','none','FontSize',12);

annotation(f,'textbox',[0.72 0.26 0.07 0.06],...
    'String',{'Forest loss\newline    pixels'},'LineStyle','none','FontSize',11);

annotation(f,'textbox',[0.78 0.26 0.07 0.06],...
    'String',{'Forest gain\newline    pixels'},'LineStyle','none','FontSize',11);

%% Deseaon & Detrend
% Loading data
load('D:\assignment\ForestResilience\response\de_season_trend\output\Tfcc\TfccTres')

IND = find((K==1)&(FC_avg>5));

FC_diff_b = nan(3600,7200);
FC_diff_b(IND) = FC_diff(IND);

edges = -0.25:0.02:0.25;
[N,edges,bin] = histcounts(FC_diff_b,edges);
X = -0.24:0.02:0.24;
Y = nan(length(edges)-1,1);
V = nan(length(edges)-1,1);
SEM = nan(length(edges)-1,1);
err = nan(length(edges)-1,2);

for i = 1:length(edges)-1
    
    ind = find(bin==i);
    Y(i) = wmean(Tfcc(ind),W(ind),'omitnan');
    V(i) = (nanvar(Tfcc(ind),W(ind))^(0.5));
    SEM(i) = V(i)/sqrt(length(Tfcc(ind))); 
    ts = tinv([0.025  0.975],length(Tfcc(ind))-1); 
    err(i,:) = ts*SEM(i); 

end

indl = find((K==1)&(FC_avg>0.05)&(FC_diff_b<0));
indg = find((K==1)&(FC_avg>0.05)&(FC_diff_b>0));

dataX1 = FC_diff_b(indl);
dataY1 = Tfcc(indl);
[xData1, yData1] = prepareCurveData( dataX1, dataY1 );
[fitresult1, gof1] = fit( xData1, yData1, 'poly1' );

dataX2 = FC_diff_b(indg);
dataY2 = Tfcc(indg);
[xData2, yData2] = prepareCurveData( dataX2, dataY2 );
[fitresult2, gof2] = fit( xData2, yData2, 'poly1' );

% figure g,h
hg = axes('Position',[0.125 0.073 0.5 0.17]);

X1 = X(1:13); X2 = X(13:end);
Y1 = Y(1:13); Y2 = Y(13:end);
b1 = bar(X1,Y1,0.8);
b1.FaceColor = mycolor(1,:);
b1.FaceAlpha = 0.5;
b1.EdgeColor = "none";
box on
hold on
b2 = bar(X2,Y2,0.8);
b2.FaceColor = mycolor(2,:);
b2.FaceAlpha = 0.5;
b2.EdgeColor = "none"; 
box on
hold on
er1 = errorbar(X1,Y1,err(1:13,1),err(1:13,2), ...
    'Linestyle','None','LineWidth',1,'CapSize',5,'Color',mycolor(1,:));
er2 = errorbar(X2,Y2,err(13:end,1),err(13:end,2), ...
    'Linestyle','None','LineWidth',1,'CapSize',5,'Color',mycolor(2,:));
ylim([-0.08, 0.08]);

plot([0,0],[-0.08,0.08],'k','Linewidth',0.5,'LineStyle','--')

ylabel({'\DeltaTAC of','kNDVI'});
xlabel({'\DeltaFC derived from MOD44B (%)'},'FontSize',14);
set(hg,'FontSize',14,'GridLineWidth',0.7,'LineWidth',1,'MinorGridAlpha',...
    0.2,'XMinorGrid','on','XTick',[-0.25 -0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25],'XTickLabel',...
    {'-25','-20','-15','-10','-5','0','5','10','15','20','25'},'YMinorGrid','on', ...
    'YTick',[-0.08 -0.04 0 0.04 0.08],'YTickLabel',{'-0.08','-0.04','0','0.04','0.08'});

hh = axes('Position',[0.707 0.073 0.16 0.17]);

x=[1,2];   
y=[fitresult1.p1,fitresult2.p1];

err1 = [ 0.012,  0.010];
b3 = bar(x(1),y(1),0.5);
hold on 
er3 = errorbar(x(1),y(1),err1(1));
b3.FaceColor = mycolor2(1,:);
b3.FaceAlpha = 0.7;
b3.EdgeColor = 'none';
b3.LineWidth = 1;
er3.LineWidth = 1.5;
er3.Color = mycolor2(1,:);
er3.CapSize = 10;

b4 = bar(x(2),y(2),0.5);
hold on 
er4 = errorbar(x(2),y(2),err1(2));
b4.FaceColor = mycolor2(2,:);
b4.FaceAlpha = 0.7;
b4.EdgeColor = 'none';
b4.LineWidth = 1;
er4.LineWidth = 1.5;
er4.Color = mycolor2(2,:);
er4.CapSize = 10;

ylim([-0.3, 0]);
set(gca, 'YDir', 'reverse');

ylabel(['Regional sensitivity    ';'of \DeltaTAC to \DeltaFC']);
set(hh,'FontSize',13,'LineWidth',1,'XTick',[1 2],...
    'XTickLabel',{'',''},'YTickLabel',{'-3','-2','-1','0'},'Box','off');

annotation(f,'textbox',[0.706 0.213 0.051 0.04],'String','(×10^-^3)',...
    'LineStyle','none','FontSize',12);

annotation(f,'textbox',[0.72 0.016 0.07 0.06],...
    'String',{'Forest loss\newline    pixels'},'LineStyle','none','FontSize',11);

annotation(f,'textbox',[0.78 0.016 0.07 0.06],...
    'String',{'Forest gain\newline    pixels'},'LineStyle','none','FontSize',11);

%% VOD
load('D:\assignment\ForestResilience\vcf\output\FC\FC_diff\FC_diff_025vod')
load('D:\assignment\ForestResilience\response\vod\output\Tfcc\output\TfccTres_new')
load('D:\assignment\ForestResilience\koppen\output\Koppen_025')
load('D:\assignment\ForestResilience\vcf\output\whole\FC_avg_025')

IND = find((K_025==1)&(FC_avg_025>5));

FC_diff_b = nan(720,1440);
FC_diff_b(IND) = (FC_diff_025vod(IND))./100;

edges = -0.10:0.005:0.10;
[N,edges,bin] = histcounts(FC_diff_b,edges);
X = -0.11:0.005:0.11;
Y = nan(length(edges)-1,1);
V = nan(length(edges)-1,1);
SEM = nan(length(edges)-1,1);
err = nan(length(edges)-1,2);

clear lat lon LON LAT
cellsize = 0.25;
lat(:,1) = -90+cellsize/2:cellsize:90-cellsize/2;
lon(1,:) = -180+cellsize/2:cellsize:180-cellsize/2;
[LON,LAT] = meshgrid(lon,flipud(lat));
W = cdtarea(LAT,LON,'km2');

for i = 1:length(edges)-1
    
    ind = find(bin==i);
    Y(i) = wmean(Tfcc(ind),W(ind),'omitnan');
    V(i) = (nanvar(Tfcc(ind),W(ind))^(0.5));
    SEM(i) = V(i)/sqrt(length(Tfcc(ind))); 
    ts = tinv([0.025  0.975],length(Tfcc(ind))-1); 
    err(i,:) = ts*SEM(i); 

end

indl = find((K_025==1)&(FC_avg_025>5)&(FC_diff_b<0));
indg = find((K_025==1)&(FC_avg_025>5)&(FC_diff_b>0));

dataX1 = FC_diff_b(indl);
dataY1 = Tfcc(indl);
[xData1, yData1] = prepareCurveData( dataX1, dataY1 );
[fitresult1, gof1] = fit( xData1, yData1, 'poly1' );

dataX2 = FC_diff_b(indg);
dataY2 = Tfcc(indg);
[xData2, yData2] = prepareCurveData( dataX2, dataY2 );
[fitresult2, gof2] = fit( xData2, yData2, 'poly1' );

% figure a,b
ha = axes('Position',[0.125,0.812,0.5,0.17]);
X1 = X(3:22); X2 = X(24:43);
Y1 = Y(1:20); Y2 = Y(21:end);
b1 = bar(X1,Y1,0.8);
b1.FaceColor = mycolor(1,:);
b1.FaceAlpha = 0.5;
b1.EdgeColor = "none";
box on
hold on
b2 = bar(X2,Y2,0.8);
b2.FaceColor = mycolor(2,:);
b2.FaceAlpha = 0.5;
b2.EdgeColor = "none";
box on
hold on
er1 = errorbar(X1,Y1,err(1:20,1),err(1:20,2), ...
    'Linestyle','None','LineWidth',1,'CapSize',5,'Color',mycolor(1,:));
er2 = errorbar(X2,Y2,err(21:end,1),err(21:end,2), ...
    'Linestyle','None','LineWidth',1,'CapSize',5,'Color',mycolor(2,:));
ylim([-0.08, 0.08]);

plot([0,0],[-0.08,0.08],'k','Linewidth',0.5,'LineStyle','--')

ylabel({'\DeltaTAC of','VOD'});
xlabel({'\DeltaFC derived from MOD44B (%)'},'FontSize',14);

set(ha,'FontSize',14,'GridLineWidth',0.7,'LineWidth',1,'MinorGridAlpha',...
    0.2,'XMinorGrid','on','XTick',...
    [-0.1 -0.08 -0.06 -0.04 -0.02 0 0.02 0.04 0.06 0.08 0.1],'XTickLabel',...
    {'-10','-8','-6','-4','-2','0','2','4','6','8','10'},'YMinorGrid','on',...
    'YTick',[-0.08 -0.04 0 0.04 0.08],'YTickLabel',...
    {'-0.08','-0.04','0','0.04','0.08'});

% Legend
annotation(f,'rectangle',[0.53 0.95 0.016 0.015],...
    'LineStyle','none','FaceColor',[0.70 0.59 0.47]);
annotation(f,'rectangle',[0.53 0.93 0.016 0.015],...
    'LineStyle','none','FaceColor',[0.50 0.68 0.658]);

annotation(f,'textbox',[0.55 0.937 0.08 0.04],...
    'String',{'Forest loss'},'LineStyle','none','FontSize',12);

annotation(f,'textbox',[0.55 0.916 0.08 0.04],...
    'String',{'Forest gain'},'LineStyle','none','FontSize',12);

hb = axes('Position',[0.707 0.812 0.16 0.17]);
x=[1,2];   
y=[fitresult1.p1,fitresult2.p1];

err1 = [ 0.076,  0.077];
b3 = bar(x(1),y(1),0.5);
hold on 
er3 = errorbar(x(1),y(1),err1(1));
b3.FaceColor = mycolor2(1,:);
b3.FaceAlpha = 0.7;
b3.EdgeColor = 'none';
b3.LineWidth = 1;
er3.LineWidth = 1.5;
er3.Color = mycolor2(1,:);
er3.CapSize = 10;

b4 = bar(x(2),y(2),0.5);
hold on 
er4 = errorbar(x(2),y(2),err1(2));
b4.FaceColor = mycolor2(2,:);
b4.FaceAlpha = 0.7;
b4.EdgeColor = 'none';
b4.LineWidth = 1;
er4.LineWidth = 1.5;
er4.Color = mycolor2(2,:);
er4.CapSize = 10;

ylim([-0.6, 0]);
set(gca, 'YDir', 'reverse');

ylabel({'Regional sensitivity ','of \DeltaTAC to \DeltaFC'},'FontSize',13);
set(hb,'FontSize',13,'LineWidth',1,'XTick',[1 2],...
    'XTickLabel',{'',''},'YTick',[-0.8 -0.6 -0.4 -0.2 0],'YTickLabel',...
    {'-8','-6','-4','-2','0'},'Box','off');

annotation(f,'textbox',[0.706 0.95 0.06 0.04],...
    'String',{'(×10^-^3)'},'LineStyle','none','FontSize',12);

annotation(f,'textbox',[0.72 0.75 0.07 0.06],...
    'String',{'Forest loss\newline    pixels'},'LineStyle','none','FontSize',11);

annotation(f,'textbox',[0.78 0.75 0.07 0.06],...
    'String',{'Forest gain\newline    pixels'},'LineStyle','none','FontSize',11);