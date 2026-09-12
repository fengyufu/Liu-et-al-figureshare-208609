%% Supplementary Fig. 4. Probability distribution of sensitivity distinguished by intact, regenerating, and planted forest.

%% 

clear
clc;

%% Plotting histograms of local sensitivity at different scales (0.25 degree)

load([pwd,'\input\data\alef_25'])
load([pwd,'\input\data\FML_25'])
load([pwd,'\input\data\FC_diff_25'])

% Creating grid
cellsize = 0.05;
lat(:,1) = -90+cellsize/2:cellsize:90-cellsize/2;
lon(1,:) = -180+cellsize/2:cellsize:180-cellsize/2;
[LON,LAT] = meshgrid(lon,flipud(lat));
W = cdtarea(LAT,LON,'km2');

ind_pri_l = find((FML_t == 11)&(isfinite(alef))&(FC_diff_t < 0));
ind_pri_g = find((FML_t == 11)&(isfinite(alef))&(FC_diff_t > 0));

ind_rege_l = find((FML_t == 20)&(isfinite(alef))&(FC_diff_t < 0));
ind_rege_g = find((FML_t == 20)&(isfinite(alef))&(FC_diff_t > 0));

ind_plan_l = find(((FML_t == 31)|(FML_t == 32))&(isfinite(alef))&(FC_diff_t < 0));
ind_plan_g = find(((FML_t == 31)|(FML_t == 32))&(isfinite(alef))&(FC_diff_t > 0));

diff1 = wmean(alef(ind_pri_l),W(ind_pri_l),'omitnan') - wmean(alef(ind_pri_g),W(ind_pri_g),'omitnan');  [h1,p1] = ttest2(alef(ind_pri_l),alef(ind_pri_g),"Alpha",0.05,"Tail","both");
diff2 = wmean(alef(ind_rege_l),W(ind_rege_l),'omitnan') - wmean(alef(ind_rege_g),W(ind_rege_g),'omitnan'); [h2,p2] = ttest2(alef(ind_rege_l),alef(ind_rege_g),"Alpha",0.05,"Tail","both");
diff3 = wmean(alef(ind_plan_l),W(ind_plan_l),'omitnan') - wmean(alef(ind_plan_g),W(ind_plan_g),'omitnan'); [h3,p3] = ttest2(alef(ind_plan_l),alef(ind_plan_g),"Alpha",0.05,"Tail","both");


f = figure('Position',[0 0 1700 840]);

ha = subplot(3,3,1);
histogram(alef(ind_pri_l),'DisplayName','Forest loss','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.80,0.47,0.75],'Normalization','pdf','NumBins',200);
hold on
histogram(alef(ind_pri_g),'DisplayName','Forest gain','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.30,0.64,0.64],'Normalization','pdf','NumBins',200);

set(ha,'FontSize',14,'GridLineWidth',1,'LineWidth',1.5,'TickDir','out','box','off',...
'xlim',[-2 2],'ylim',[0 1], ...
'XTickLabel',{'-20','-10','0','10','20'},'YTickLabel',{'0','0.5','1.0'});
title({'IF       -0.0010(*)'});


hb = subplot(3,3,2);
histogram(alef(ind_rege_l),'DisplayName','Forest loss','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.80,0.47,0.75],'Normalization','pdf','NumBins',200);
hold on
histogram(alef(ind_rege_g),'DisplayName','Forest gain','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.30,0.64,0.64],'Normalization','pdf','NumBins',200);

set(hb,'FontSize',14,'GridLineWidth',1,'LineWidth',1.5,'TickDir','out','box','off',...
'xlim',[-2 2],'ylim',[0 1], ...
'XTickLabel',{'-20','-10','0','10','20'},'YTickLabel',{'0','0.5','1.0'});
title({'RF       -0.0003(*)'});


hc = subplot(3,3,3);
histogram(alef(ind_plan_l),'DisplayName','Forest loss','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.80,0.47,0.75],'Normalization','pdf','NumBins',250);
hold on
histogram(alef(ind_plan_g),'DisplayName','Forest gain','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.30,0.64,0.64],'Normalization','pdf','NumBins',250);

set(hc,'FontSize',14,'GridLineWidth',1,'LineWidth',1.5,'TickDir','out','box','off',...
'xlim',[-2 2],'ylim',[0 1], ...
'XTickLabel',{'-20','-10','0','10','20'},'YTickLabel',{'0','0.5','1.0'});
title({'PF       -0.0007(*)'});


annotation(f,'textbox',[0.91 0.70 0.04 0.04],'String',{'(x10^-^3)'},...
    'LineStyle','none','FontSize',12); % convertion into sensitivty of ΔTAC to 1% forest cover change (x10^-2)

%% Plotting histograms of local sensitivity at different scales (0.5 degree)

load([pwd,'\input\data\alef_50'])
load([pwd,'\input\data\FML_50'])
load([pwd,'\input\data\FC_diff_50'])

ind_pri_l = find((FML_t == 11)&(isfinite(alef))&(FC_diff_t < 0));
ind_pri_g = find((FML_t == 11)&(isfinite(alef))&(FC_diff_t > 0));

ind_rege_l = find((FML_t == 20)&(isfinite(alef))&(FC_diff_t < 0));
ind_rege_g = find((FML_t == 20)&(isfinite(alef))&(FC_diff_t > 0));

ind_plan_l = find(((FML_t == 31)|(FML_t == 32))&(isfinite(alef))&(FC_diff_t < 0));
ind_plan_g = find(((FML_t == 31)|(FML_t == 32))&(isfinite(alef))&(FC_diff_t > 0));

diff1 = wmean(alef(ind_pri_l),W(ind_pri_l),'omitnan') - wmean(alef(ind_pri_g),W(ind_pri_g),'omitnan');  [h1,p1] = ttest2(alef(ind_pri_l),alef(ind_pri_g),"Alpha",0.05,"Tail","both");
diff2 = wmean(alef(ind_rege_l),W(ind_rege_l),'omitnan') - wmean(alef(ind_rege_g),W(ind_rege_g),'omitnan'); [h2,p2] = ttest2(alef(ind_rege_l),alef(ind_rege_g),"Alpha",0.05,"Tail","both");
diff3 = wmean(alef(ind_plan_l),W(ind_plan_l),'omitnan') - wmean(alef(ind_plan_g),W(ind_plan_g),'omitnan'); [h3,p3] = ttest2(alef(ind_plan_l),alef(ind_plan_g),"Alpha",0.05,"Tail","both");

hd = subplot(3,3,4);
histogram(alef(ind_pri_l),'DisplayName','Forest loss','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.80,0.47,0.75],'Normalization','pdf','NumBins',150);
hold on
histogram(alef(ind_pri_g),'DisplayName','Forest gain','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.30,0.64,0.64],'Normalization','pdf','NumBins',200);

set(hd,'FontSize',14,'GridLineWidth',1,'LineWidth',1.5,'TickDir','out','box','off',...
'xlim',[-1.5 1.5],'ylim',[0 1.5], ...
'XTickLabel',{'-15','-10','-5','0','5','10','15'},'YTickLabel',{'0','0.5','1.0','1.5'});
title({'IF       -0.0007(*)'});
ylabel({'Probability density'},'FontSize',16);


he = subplot(3,3,5);
histogram(alef(ind_rege_l),'DisplayName','Forest loss','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.80,0.47,0.75],'Normalization','pdf','NumBins',200);
hold on
histogram(alef(ind_rege_g),'DisplayName','Forest gain','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.30,0.64,0.64],'Normalization','pdf','NumBins',250);

set(he,'FontSize',14,'GridLineWidth',1,'LineWidth',1.5,'TickDir','out','box','off',...
'xlim',[-1.5 1.5],'ylim',[0 1.5], ...
'XTickLabel',{'-15','-10','-5','0','5','10','15'},'YTickLabel',{'0','0.5','1.0','1.5'});
title({'RF       -0.0000(*)'});

hf = subplot(3,3,6);
histogram(alef(ind_plan_l),'DisplayName','Forest loss','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.80,0.47,0.75],'Normalization','pdf','NumBins',200);
hold on
histogram(alef(ind_plan_g),'DisplayName','Forest gain','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.30,0.64,0.64],'Normalization','pdf','NumBins',200);

set(hf,'FontSize',14,'GridLineWidth',1,'LineWidth',1.5,'TickDir','out','box','off',...
'xlim',[-1.5 1.5],'ylim',[0 1.5], ...
'XTickLabel',{'-15','-10','-5','0','5','10','15'},'YTickLabel',{'0','0.5','1.0','1.5'});
title({'PF       -0.0008(*)'});

annotation(f,'textbox',[0.91 0.40 0.04 0.04],'String',{'(x10^-^3)'},...
    'LineStyle','none','FontSize',12); % convertion into sensitivty of ΔTAC to 1% forest cover change (x10^-2)

%% Plotting histograms of local sensitivity at different scales (2 degree)

load([pwd,'\input\data\alef_200'])
load([pwd,'\input\data\FML_200'])
load([pwd,'\input\data\FC_diff_200'])

ind_pri_l = find((FML_t == 11)&(isfinite(alef))&(FC_diff_t < 0));
ind_pri_g = find((FML_t == 11)&(isfinite(alef))&(FC_diff_t > 0));

ind_rege_l = find((FML_t == 20)&(isfinite(alef))&(FC_diff_t < 0));
ind_rege_g = find((FML_t == 20)&(isfinite(alef))&(FC_diff_t > 0));

ind_plan_l = find(((FML_t == 31)|(FML_t == 32))&(isfinite(alef))&(FC_diff_t < 0));
ind_plan_g = find(((FML_t == 31)|(FML_t == 32))&(isfinite(alef))&(FC_diff_t > 0));

diff1 = wmean(alef(ind_pri_l),W(ind_pri_l),'omitnan') - wmean(alef(ind_pri_g),W(ind_pri_g),'omitnan');  [h1,p1] = ttest2(alef(ind_pri_l),alef(ind_pri_g),"Alpha",0.05,"Tail","both");
diff2 = wmean(alef(ind_rege_l),W(ind_rege_l),'omitnan') - wmean(alef(ind_rege_g),W(ind_rege_g),'omitnan'); [h2,p2] = ttest2(alef(ind_rege_l),alef(ind_rege_g),"Alpha",0.05,"Tail","both");
diff3 = wmean(alef(ind_plan_l),W(ind_plan_l),'omitnan') - wmean(alef(ind_plan_g),W(ind_plan_g),'omitnan'); [h3,p3] = ttest2(alef(ind_plan_l),alef(ind_plan_g),"Alpha",0.05,"Tail","both");


hg = subplot(3,3,7);
histogram(alef(ind_pri_l),'DisplayName','Forest loss','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.80,0.47,0.75],'Normalization','pdf','NumBins',100);
hold on
histogram(alef(ind_pri_g),'DisplayName','Forest gain','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.30,0.64,0.64],'Normalization','pdf','NumBins',250);

set(hg,'FontSize',14,'GridLineWidth',1,'LineWidth',1.5,'TickDir','out','box','off',...
'xlim',[-1 1],'ylim',[0 2.5], ...
'XTickLabel',{'-10','-5','0','5','10'},'YTickLabel',{'0','0.5','1.0','1.5','2.0','2.5'});
title({'IF       -0.0009(*)'});


hh = subplot(3,3,8);
histogram(alef(ind_rege_l),'DisplayName','Forest loss','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.80,0.47,0.75],'Normalization','pdf','NumBins',80);
hold on
histogram(alef(ind_rege_g),'DisplayName','Forest gain','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.30,0.64,0.64],'Normalization','pdf','NumBins',100);

set(hh,'FontSize',14,'GridLineWidth',1,'LineWidth',1.5,'TickDir','out','box','off',...
'xlim',[-1 1],'ylim',[0 2.5], ...
'XTickLabel',{'-10','-5','0','5','10'},'YTickLabel',{'0','0.5','1.0','1.5','2.0','2.5'});
title({'RF        0.0002(*)'});
xlabel({'Sensitivity of \DeltaTAC to \DeltaFC'},'FontSize',16);


hi = subplot(3,3,9);
histogram(alef(ind_plan_l),'DisplayName','Forest loss','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.80,0.47,0.75],'Normalization','pdf','NumBins',50);
hold on
histogram(alef(ind_plan_g),'DisplayName','Forest gain','FaceAlpha',0.55,...
    'EdgeColor',[1 1 1],'FaceColor',[0.30,0.64,0.64],'Normalization','pdf','NumBins',100);

set(hi,'FontSize',14,'GridLineWidth',1,'LineWidth',1.5,'TickDir','out','box','off',...
'xlim',[-1 1],'ylim',[0 2.5], ...
'XTickLabel',{'-10','-5','0','5','10'},'YTickLabel',{'0','0.5','1.0','1.5','2.0','2.5'});
title({'PF        0.0015(*)'});

annotation(f,'textbox',[0.91 0.10 0.04 0.04],'String',{'(x10^-^3)'},...
    'LineStyle','none','FontSize',12); % convertion into sensitivty of ΔTAC to 1% forest cover change (x10^-2)
