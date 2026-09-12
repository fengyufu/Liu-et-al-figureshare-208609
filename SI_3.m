%% Supplementary Fig. 3. GPP and the cover fraction of local sensitivity under different forest types.

%% 

clear
clc;

%% Loading data

load([pwd,'\input\data\Koppen_BF'])
load([pwd,'\input\data\GPP_multiyear_avg'])
load([pwd,'\input\data\FML'])
load([pwd,'\input\data\alef_100']) % local sensitivity calculated in 1d x 1d window, 
                                   % TAC change per 100% forest cover change

%% Figure plotting 

f = figure('Position',[0 0 900 500]);

% Creating grid
cellsize = 0.05;
lat(:,1) = -90+cellsize/2:cellsize:90-cellsize/2;
lon(1,:) = -180+cellsize/2:cellsize:180-cellsize/2;
[LON,LAT] = meshgrid(lon,flipud(lat));
W = cdtarea(LAT,LON,'km2');

ind_pri_005 = find(FML == 11);
ind_rege_005 = find(FML == 20);
ind_plan_005 = find((FML == 31)|(FML == 32));

data = nan(3,1); 
data(1,1) = wmean(GPP_avg(ind_pri_005),W(ind_pri_005),'omitnan');
data(2,1) = wmean(GPP_avg(ind_rege_005),W(ind_rege_005),'omitnan');
data(3,1) = wmean(GPP_avg(ind_plan_005),W(ind_plan_005),'omitnan');

var = nan(3,1);
var(1,1) = std(GPP_avg(ind_pri_005),W(ind_pri_005),'omitnan');
var(2,1) = std(GPP_avg(ind_rege_005),W(ind_rege_005),'omitnan');
var(3,1) = std(GPP_avg(ind_plan_005),W(ind_plan_005),'omitnan');

% Figure a
ha = axes('Position',[0.133 0.312 0.296 0.522]);

b1 = bar([1 2 3],data,'FaceAlpha',1,'FaceColor',[0.76,0.82,0.57],...
    'EdgeColor',[0.76,0.81,0.64],'BarWidth',0.6);
hold on
er1 = errorbar([1 2 3],data,var,'LineStyle','none','LineWidth',1.5,...
    'Color',[0.42,0.52,0.13],'CapSize',12);
ylim([0, 1.5]);
set(ha,'FontSize',12,'LineWidth',1,'TickDir','out','XTick',[1 2 3],...
    'XTickLabel',{'IF','RF','PF'});
ylabel({'Mean annual GPP','(kgC m^-² year^-^1)'});


% Figure b

hb = axes('Position',[0.526 0.305 0.294 0.521]);

load([pwd,'\input\data\FML_1d'])
ind_pri = find(FML_t == 11);
ind_rege = find(FML_t == 20);
ind_plan = find((FML_t == 31)|(FML_t == 32));

pro_pos_pri = length(find((alef(ind_pri)>0)&(p(ind_pri)<0.05)))./(length(find((alef(ind_pri)<0)&(p(ind_pri)<0.05)))+length(find((alef(ind_pri)>0)&(p(ind_pri)<0.05))));
pro_neg_pri = 1-pro_pos_pri;

pro_pos_rege = length(find((alef(ind_rege)>0)&(p(ind_rege)<0.05)))./(length(find((alef(ind_rege)<0)&(p(ind_rege)<0.05)))+length(find((alef(ind_rege)>0)&(p(ind_rege)<0.05))));
pro_neg_rege = 1-pro_pos_rege;

pro_pos_plan = length(find((alef(ind_plan)>0)&(p(ind_plan)<0.05)))./(length(find((alef(ind_plan)<0)&(p(ind_plan)<0.05)))+length(find((alef(ind_plan)>0)&(p(ind_plan)<0.05))));
pro_neg_plan = 1-pro_pos_plan;

% pro_pos_pri = length(find((alef(ind_pri)>0)))./(length(find(isfinite(alef(ind_pri)))));
% pro_neg_pri = 1-pro_pos_pri;
% 
% pro_pos_rege = length(find((alef(ind_rege)>0)))./(length(find(isfinite(alef(ind_rege)))));
% pro_neg_rege = 1-pro_pos_rege;
% 
% pro_pos_plan = length(find((alef(ind_plan)>0)))./(length(find(isfinite(alef(ind_plan)))));
% pro_neg_plan = 1-pro_pos_plan;

y = [pro_pos_pri pro_neg_pri; pro_pos_rege pro_neg_rege; pro_pos_plan pro_neg_plan];
b = bar(y,'BarWidth',0.7,'BarLayout','stacked');

set(b(1),'DisplayName','sensitivity < 0','FaceAlpha', 0.75,'FaceColor',[0.82 0.33 0.27]);
set(b(2),'DisplayName','sensitivity > 0','FaceColor',[0.51 0.73 0.85]);

ylabel({'Cover fraction'},'FontSize',14);

box(hb,'on');

set(hb,'FontSize',12,'LineWidth',1,'TickDir','out','XTick',[1 2 3],...
    'XTickLabel',{'IF','RF','PF'},'YTick',[0 0.2 0.4 0.6 0.8 1]);
lgd= legend(hb,'show');
set(lgd,'Position',[0.83 0.31 0.16 0.09],'FontSize',12,'EdgeColor',[1 1 1]);
