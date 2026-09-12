%% Supplementary Fig. 2. Analysis of forest resilience and its sensitivity to forest cover change across different forest types.

%% 

clear
clc;

%% Loading data

load([pwd,'\input\data\long-term_TAC'])
load([pwd,'\input\data\FC_multiyear_avg'])
load([pwd,'\input\data\FC_diff_3y'])
load([pwd,'\input\data\Koppen_BF'])
load([pwd,'\input\data\Tfcc_BF'],'Tfcc')
load([pwd,'\input\data\FML'])

%% Selecting pixels in boreal region

indnan = find((K~=1)|(FC_avg <5));
ACF(indnan) = nan;
FC_avg(indnan) = nan;
Tfcc(indnan) = nan;

%% Categorizing FC and long-term TAC into bins

xlimset = [5, 60; 0, 0.7];
xx = FC_avg;
yy = ACF;

minX = xlimset(1,1);
maxX = xlimset(1,2);
minY = xlimset(2,1);
maxY = xlimset(2,2);

xx(xx<minX) = nan;
xx(xx>maxX) = nan;
yy(yy<minY) = nan;
yy(yy>maxY) = nan;

[INDnan] = find(isnan(xx)|isnan(yy));
xx(INDnan) = [];
yy(INDnan) = [];

targetSize = [60 70];

xxBin = round( (xx-minX)/(maxX-minX)*(targetSize(1)-1) ) +1;
yyBin = round( (yy-minY)/(maxY-minY)*(targetSize(2)-1) ) +1;
Z = nan(targetSize(1), targetSize(2));

for ibin = 1:targetSize(1)
    for jbin=1:targetSize(2)
        
        INDbin = find((xxBin==ibin)&(yyBin==jbin));
        Z(ibin, jbin) = length(INDbin);
        
    end
end

Z = (Z./length(xx)).*100; 

%% Computing long-term TAC in each FC bin across different forest cover types

% Creating grid
cellsize = 0.05;
lat(:,1) = -90+cellsize/2:cellsize:90-cellsize/2;
lon(1,:) = -180+cellsize/2:cellsize:180-cellsize/2;
[LON,LAT] = meshgrid(lon,flipud(lat));
W = cdtarea(LAT,LON,'km2');

ind_pri_005 = find((FML == 11)&(K==1)&(FC_avg > 5));
ind_rege_005 = find((FML == 20)&(K==1)&(FC_avg > 5));
ind_plan_005 = find((FML == 31)&(K==1)&(FC_avg > 5));

pri_b = nan(3600,7200);
pri_b(ind_pri_005) = ACF((ind_pri_005));

rege_b = nan(3600,7200);
rege_b(ind_rege_005) = ACF((ind_rege_005));

plan_b = nan(3600,7200);
plan_b(ind_plan_005) = ACF((ind_plan_005));


IND = find((K==1)&(FC_avg>5));

FC_avg_b = nan(3600,7200);
FC_avg_b(IND) = FC_avg(IND);

edges = 5:1:60;
[N,edges,bin] = histcounts(FC_avg_b,edges);
X = 5:1:60;

Y_pri = nan(length(edges)-1,1);
for i = 1:length(edges)-1
    
    ind = find(bin==i);
    Y_pri(i) = wmean(pri_b(ind),W(ind),'omitnan');
  
end
Y_pri = movmean(Y_pri,3);


Y_rege = nan(length(edges)-1,1);
for i = 1:length(edges)-1
    
    ind = find(bin==i);
    Y_rege(i) = wmean(rege_b(ind),W(ind),'omitnan');
  
end
Y_rege = movmean(Y_rege,3);


Y_plan = nan(length(edges)-1,1);
for i = 1:length(edges)-1
    
    ind = find(bin==i);
    Y_plan(i) = wmean(plan_b(ind),W(ind),'omitnan');
  
end
Y_plan = movmean(Y_plan,3);

%% Figure plotting (figure a)

f = figure('Position',[0 0 1000 750]);

ha = axes('Position',[0.087 0.575 0.289 0.385]);

% Re-gridding
Xvet = minX:((maxX-minX)/(targetSize(1)-1)):maxX;
Yvet = minY:((maxY-minY)/(targetSize(2)-1)):maxY;
[Xmat, Ymat] = meshgrid(Xvet, Yvet);

cmap = brewermap(12,'YlOrBr');
cmap(7:end,:) = [];

pcolor(Xmat,Ymat,Z'); shading flat; axis square;
axis([min(Xmat(:)) max(Xmat(:)) min(Ymat(:)) max(Ymat(:))])
set(gca,'layer','top')
colormap(gca,cmap)
clim([0.01,0.07])
colorbar(ha,"Position",[0.395,0.571,0.021,0.385])
set(ha,'FontSize',12,'Layer','top','LineWidth',1)
hold on 

plot(X(1:end-1),Y_pri,'LineWidth',2.2,'Color',[0.37,0.65,0.61])
plot(X(1:end-1),Y_rege,'LineWidth',2.2,'Color',[0.67,0.81,0.36])
plot(X(1:end-1),Y_plan,'LineWidth',2.2,'Color',[0.64,0.47,0.62])

ylabel({'Long-term TAC'},'FontSize',13);
xlabel({'FC (%)'});

%% Computing regional sensitivity across different forest types for both loss and gain

indnan = find((K==0)|(FML == 53));
FML(indnan) = nan;

ind_pri_l_005 = find((FML == 11)&(K==1)&(FC_avg > 5)&(FC_diff < 0));
ind_pri_g_005 = find((FML == 11)&(K==1)&(FC_avg > 5)&(FC_diff > 0));
ind_rege_l_005 = find((FML == 20)&(K==1)&(FC_avg > 5)&(FC_diff < 0));
ind_rege_g_005 = find((FML == 20)&(K==1)&(FC_avg > 5)&(FC_diff > 0));
ind_plan_l_005 = find(((FML == 31)|(FML == 32))&(K==1)&(FC_avg > 5)&(FC_diff < 0));
ind_plan_g_005 = find(((FML == 31)|(FML == 32))&(K==1)&(FC_avg > 5)&(FC_diff > 0));

dataX1_l = FC_diff(ind_pri_l_005);
dataY1_l = Tfcc(ind_pri_l_005);
[xData, yData] = prepareCurveData( dataX1_l, dataY1_l );
[fitresult1_l, gof1_l] = fit( xData, yData, 'poly1' );
mdl1_l = fitlm(dataX1_l, dataY1_l, 'linear' );

dataX1_g = FC_diff(ind_pri_g_005);
dataY1_g = Tfcc(ind_pri_g_005);
[xData, yData] = prepareCurveData( dataX1_g, dataY1_g );
[fitresult1_g, gof1_g] = fit( xData, yData, 'poly1' );
mdl1_g = fitlm( dataX1_g, dataY1_g, 'linear' );

dataX2_l = FC_diff(ind_rege_l_005);
dataY2_l = Tfcc(ind_rege_l_005);
[xData, yData] = prepareCurveData( dataX2_l, dataY2_l );
[fitresult2_l, gof2_l] = fit( xData, yData, 'poly1' );
mdl2_l = fitlm( dataX2_l, dataY2_l, 'linear' );

dataX2_g = FC_diff(ind_rege_g_005);
dataY2_g = Tfcc(ind_rege_g_005);
[xData, yData] = prepareCurveData( dataX2_g, dataY2_g );
[fitresult2_g, gof2_g] = fit( xData, yData, 'poly1' );
mdl2_g = fitlm( dataX2_g, dataY2_g, 'linear' );

dataX3_l = FC_diff(ind_plan_l_005);
dataY3_l = Tfcc(ind_plan_l_005);
[xData, yData] = prepareCurveData( dataX3_l, dataY3_l );
[fitresult3_l, gof3_l] = fit( xData, yData, 'poly1' );
mdl3_l = fitlm( dataX3_l, dataY3_l, 'linear' );

dataX3_g = FC_diff(ind_plan_g_005);
dataY3_g = Tfcc(ind_plan_g_005);
[xData, yData] = prepareCurveData( dataX3_g, dataY3_g );
[fitresult3_g, gof3_g] = fit( xData, yData, 'poly1' );
mdl3_g = fitlm( dataX3_g, dataY3_g, 'linear' );

yfit1_l = fitresult1_l.p1.*dataX1_l;
yfit1_g = fitresult1_g.p1.*dataX1_g;

yfit2_l = fitresult2_l.p1.*dataX2_l;
yfit2_g = fitresult2_g.p1.*dataX2_g;

yfit3_l = fitresult3_l.p1.*dataX3_l;
yfit3_g = fitresult3_g.p1.*dataX3_g;

%% Figure plotting (figure b)

hb = axes('Position',[0.566 0.568 0.289 0.391]);

scatter(FC_diff(ind_pri_005),Tfcc(ind_pri_005),2.5,[0.7,0.75,0.75],"filled")
set(hb,'FontSize',12,'LineWidth',1,'TickLength',[0.025 0.025], 'xlim',[-80,80],'ylim',[-1,1],...
    'XTick',[-80 -40 0 40 80],'XTickLabel',{'-80','-40','0','40','80'}, ...
    'YTick',[-1 -0.5 0 0.5 1],'YTickLabel',{'-1.0','-0.5','0','0.5','1.0'});
hold on
plot(dataX1_l,yfit1_l,'Color',[0.85,0.33,0.10],'LineWidth',1.5)
plot(dataX1_g,yfit1_g,'Color',[0.85,0.33,0.10],'LineWidth',1.5)

plot([0 0],[-1 1],'LineStyle','--','Color','k');
plot([-80 80],[0 0],'LineStyle','--','Color','k');

ylabel({'\DeltaTAC'});
xlabel({'\DeltaFC (%)'});
title({'IF'},'FontWeight','bold','FontSize',14);

betal_hat = mdl1_l.Coefficients.Estimate(2);
se_beta1 = mdl1_l.Coefficients.SE(2);

betag_hat = mdl1_g.Coefficients.Estimate(2);
se_beta2 = mdl1_g.Coefficients.SE(2);

W = (betal_hat - betag_hat)^2 / (se_beta1^2 + se_beta2^2);
df = 1;
pValue = 1 - chi2cdf(W, df);

annotation(f,'textbox',[0.723 0.960 0.047 0.041],'String',{'(**)'},...
    'LineStyle','none','FontSize',12);

annotation(f,'textbox',[0.721 0.861 0.146 0.1],...
    'String',{'k_l_o_s_s= -6.0×10^-^3','k_g_a_i_n= -1.4×10^-^3'},...
    'LineStyle','none','FontSize',12);

%% Figure plotting (figure c)

hc = axes('Position',[0.088 0.081 0.289 0.380]);

scatter(FC_diff(ind_rege_005),Tfcc(ind_rege_005),2.5,[0.7,0.75,0.75],"filled")
set(hc,'FontSize',12,'LineWidth',1,'TickLength',[0.025 0.025], 'xlim',[-80,80],'ylim',[-1,1],...
    'XTick',[-80 -40 0 40 80],'XTickLabel',{'-80','-40','0','40','80'}, ...
    'YTick',[-1 -0.5 0 0.5 1],'YTickLabel',{'-1.0','-0.5','0','0.5','1.0'});

hold on
plot(dataX2_l,yfit2_l,'Color',[0.85,0.33,0.10],'LineWidth',1.5)
plot(dataX2_g,yfit2_g,'Color',[0.85,0.33,0.10],'LineWidth',1.5)

plot([0 0],[-1 1],'LineStyle','--','Color','k');
plot([-80 80],[0 0],'LineStyle','--','Color','k');

ylabel({'\DeltaTAC'});
xlabel({'\DeltaFC (%)'});
title({'RF'},'FontWeight','bold','FontSize',14);

betal_hat = mdl2_l.Coefficients.Estimate(2);
se_beta1 = mdl2_l.Coefficients.SE(2);

betag_hat = mdl2_g.Coefficients.Estimate(2);
se_beta2 = mdl2_g.Coefficients.SE(2);

W = (betal_hat - betag_hat)^2 / (se_beta1^2 + se_beta2^2);
df = 1;
pValue = 1 - chi2cdf(W, df);

annotation(f,'textbox',[0.241 0.365 0.146 0.1],...
    'String',{'k_l_o_s_s= -3.1×10^-^3','k_g_a_i_n= -1.4×10^-^3'},...
    'LineStyle','none','FontSize',12);

annotation(f,'textbox',[0.245 0.462 0.047 0.041],'String',{'(**)'},...
    'LineStyle','none','FontSize',12);

%% Figure plotting (figure d)

hd = axes('Position',[0.567 0.081 0.289 0.380]);

scatter(FC_diff(ind_plan_005),Tfcc(ind_plan_005),2.5,[0.7,0.75,0.75],"filled")
set(hd,'FontSize',12,'LineWidth',1,'TickLength',[0.025 0.025], 'xlim',[-80,80],'ylim',[-1,1],...
    'XTick',[-80 -40 0 40 80],'XTickLabel',{'-80','-40','0','40','80'}, ...
    'YTick',[-1 -0.5 0 0.5 1],'YTickLabel',{'-1.0','-0.5','0','0.5','1.0'});

hold on
plot(dataX3_l,yfit3_l,'Color',[0.85,0.33,0.10],'LineWidth',1.5)
plot(dataX3_g,yfit3_g,'Color',[0.85,0.33,0.10],'LineWidth',1.5)

plot([0 0],[-1 1],'LineStyle','--','Color','k');
plot([-80 80],[0 0],'LineStyle','--','Color','k');

ylabel({'\DeltaTAC'});
xlabel({'\DeltaFC (%)'});
title({'PF'},'FontWeight','bold','FontSize',14);

betal_hat = mdl3_l.Coefficients.Estimate(2);
se_beta1 = mdl3_l.Coefficients.SE(2);

betag_hat = mdl3_g.Coefficients.Estimate(2);
se_beta2 = mdl3_g.Coefficients.SE(2);

W = (betal_hat - betag_hat)^2 / (se_beta1^2 + se_beta2^2);
df = 1;
pValue = 1 - chi2cdf(W, df);

annotation(f,'textbox',[0.724 0.364 0.146 0.1],'String',{'k_l_o_s_s= -1.7×10^-^3','k_g_a_i_n= -1.1×10^-^3'},...
    'LineStyle','none','FontSize',12);