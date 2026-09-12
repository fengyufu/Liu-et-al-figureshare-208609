%% Supplementary Fig. 1. Relationship between forest cover and TAC in boreal forests.

%% 

clear
clc;

%% Loading data 

load([pwd,'\input\data\Koppen_BF'])
load([pwd,'\input\data\FC_multiyear_avg'])
load([pwd,'\input\data\FC_diff_3y'])
load([pwd,'\input\data\Tfcc_BF'])
load([pwd,'\input\data\long-term_TAC'])

%% Computing ΔFC and ΔTAC in each FC bin

cellsize = 0.05;
lat(:,1) = -90+cellsize/2:cellsize:90-cellsize/2;
lon(1,:) = -180+cellsize/2:cellsize:180-cellsize/2;
[LON,LAT] = meshgrid(lon,flipud(lat));
W = cdtarea(LAT,LON,'km2');


IND = find((K==1)&(FC_avg > 5));

FC_avg_b = nan(3600,7200);
FC_avg_b(IND) = FC_avg(IND);

edges = 5:1:70;
[N,edges,bin] = histcounts(FC_avg_b,edges);
X = 5:1:70;
Y1 = nan(length(edges)-1,1);
Y2 = nan(length(edges)-1,1);
V1 = nan(length(edges)-1,1);
V2 = nan(length(edges)-1,1);
SEM1 = nan(length(edges)-1,1);
SEM2 = nan(length(edges)-1,1);
Y_ci1 = nan(length(edges)-1,2);
Y_ci2 = nan(length(edges)-1,2);
err1 = nan(length(edges)-1,2);
err2 = nan(length(edges)-1,2);

for i = 1:length(edges)-1
    
    ind = find(bin==i);
    Y1(i) = wmean(FC_diff(ind),W(ind),'omitnan');
    V1(i) = (nanvar(FC_diff(ind),W(ind))^(0.5));
    SEM1(i) = V1(i)/sqrt(length(FC_diff(ind))); 
    ts1 = tinv([0.025  0.975],length(FC_diff(ind))-1); 
    err1(i,:) = ts1*SEM1(i); % Y_ci(i,:) = Y(i) + ts*SEM(i);

end

for i = 1:length(edges)-1
    
    ind = find(bin==i);
    Y2(i) = wmean(Tfcc(ind),W(ind),'omitnan');
    V2(i) = (nanvar(Tfcc(ind),W(ind))^(0.5));
    SEM2(i) = V2(i)/sqrt(length(Tfcc(ind))); 
    ts2 = tinv([0.025  0.975],length(Tfcc(ind))-1); 
    err2(i,:) = ts2*SEM2(i); % Y_ci(i,:) = Y(i) + ts*SEM(i);

end

%% Figure plotting (figure a)

f = figure('Position',[0 0 700 600]);

ha = axes('Position',[0.111 0.553 0.782 0.428]); 

yyaxis left
errorbar(X(1:end-1),Y1,err1,'LineStyle','none','Marker','o','CapSize',0,'Color',[0.51,0.84,0.95],'LineWidth',1,'MarkerSize',7) 
hold on
set(ha,'YColor',[0.30,0.75,0.93],'ylim',[-6.5,6.5],'xlim',[0,70], ...
    'YTick',[-6 -4 -2 0 2 4 6], 'YTickLabel',{'-6','-4','-2','0','2','4','6'}, ...
    'FontSize',12,'LineWidth',1);
ylabel({'\DeltaFC (%)'});

yyaxis right
errorbar(X(1:end-1),Y2,err2,'LineStyle','none','Marker','o','CapSize',0,'Color',[0.95,0.63,0.55],'LineWidth',1,'MarkerSize',7)
set(ha,'YColor',[0.90,0.50,0.33],'ylim',[-0.05,0.05],'xlim',[0,70], ...
    'YTick',[-0.05 -0.025 0 0.025 0.05], 'YTickLabel',{'-0.050','-0.025','0','0.025','0.050'});
ylabel({'\DeltaTAC'}); 

%% Figure plotting (figure b)

ind_all = find((K==1)&(FC_avg > 5));
FC_avg_all = nan(3600,7200);
FC_avg_all(ind_all) = FC_avg(ind_all);

ind_l = find((K==1)&(FC_avg > 5)&(FC_diff < 0));
FC_avg_l = nan(3600,7200);
FC_avg_l(ind_l) = FC_avg(ind_l);

ind_g = find((K==1)&(FC_avg > 5)&(FC_diff > 0));
FC_avg_g = nan(3600,7200);
FC_avg_g(ind_g) = FC_avg(ind_g);

%% 

hb = axes('Position',[0.111 0.092 0.782 0.405]); 

edges = 5:1:70;
[~,edges,bin] = histcounts(FC_avg_all,edges);

Y_all = nan(length(edges)-1,1);

for i = 1:length(edges)-1

    ind = find(bin==i);
   
    Y_all(i) = wmean(ACF(ind),W(ind),'omitnan');
    
end

Y_all = movmean(Y_all,3);

[~,edges,bin] = histcounts(FC_avg_l,edges);

Y_l = nan(length(edges)-1,1);

for i = 1:length(edges)-1

    ind = find(bin==i);
   
    Y_l(i) = wmean(ACF(ind),W(ind),'omitnan');
    
end
Y_l = movmean(Y_l,3);

[~,edges,bin] = histcounts(FC_avg_g,edges);

Y_g = nan(length(edges)-1,1);

for i = 1:length(edges)-1

    ind = find(bin==i);
   
    Y_g(i) = wmean(ACF(ind),W(ind),'omitnan');
    
end
Y_g = movmean(Y_g,3);

p = plot(X(1:end-1),Y_all,X(1:end-1),Y_l,X(1:end-1),Y_g);

set(p(1),'DisplayName','All pixels','Color',[0.651 0.651 0.651],'LineWidth',1.5);
set(p(2),'DisplayName','Forest loss pixels','Color',[0.514 0.431 0.318],'LineWidth',1.5);
set(p(3),'DisplayName','Forest gain pixels','Color',[0.325 0.655 0.647],'LineWidth',1.5);

set(hb,'FontSize',12,'LineWidth',1,'XTickLabel',{'0','10','20','30','40','50','60','70'},'YGrid','on','YTickLabel',...
    {'0.1','0.2','0.3','0.4','0.5','0.6'});

lgnd = legend(hb,'show');
set(lgnd,'EdgeColor',[1 1 1]);

ylabel({'Long-term TAC'});
xlabel({'FC (%)'});
