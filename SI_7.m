%% Supplementary Fig. 7. Land-cover type change condition of target pixels and its potential influences on regional sensitivity.

%% 

clear
clc;

%% Loading data

load([pwd,'\input\data\target_reference_pixels'])
load([pwd,'\input\data\LCT'])
load([pwd,'\input\data\Tfcc_BF'])
load([pwd,'\input\data\FC_diff_3y'])
FC_diff = FC_diff./100; % convert the range of forest cover fraction from 0-100(%) to 0-1 
load([pwd,'\input\data\Koppen_BF'])
load([pwd,'\input\data\FC_multiyear_avg'])
FC_avg = FC_avg./100;

%% Converting the projection of data for figure plotting

MASK_nochange = (LCT_2001 == LCT_2020);
ind_nochange = find((MASKt == 1)&(MASK_nochange == 1));
ind_change = find((MASKt == 1)&(MASK_nochange == 0));

NumMat = nan(3600,7200);

NumMat(ind_nochange) = 1;
NumMat(ind_change) = 2;

latlim = [40 90];
lonlim = [-180 180];

R = georasterref('RasterSize', [1000 7200], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [40 90], 'LongitudeLimits', [-180 180]);
[B,RB] = geocrop(NumMat(1:1000,:),R,latlim,lonlim);
[lat,lon] = geographicGrid(RB);

%% Recording land cover type transition

LCT_2001_sel = LCT_2001(ind_change);
LCT_2020_sel = LCT_2020(ind_change);

changeTab = zeros(14,14);
lct = [1 3 4 5 8 9 10 12 14];

for i = 1:length(lct)
    for j = 1:length(lct)
        
        ind = LCT_2001_sel == lct(i);
        
        changeTab(lct(i),lct(j)) = length(find(LCT_2020_sel(ind) == lct(j)));

    end
end


links = {'ENF', 'dnf', changeTab(1,3); 'ENF', 'dbf', changeTab(1,4); 'ENF', 'mf', changeTab(1,5); 'ENF', 'woody savannas', changeTab(1,8); 'ENF', 'savannas', changeTab(1,9); 'ENF', 'grassland', changeTab(1,10);
         'DNF', 'dbf', changeTab(3,4); 'DNF', 'mf', changeTab(3,5); 'DNF', 'woody savannas', changeTab(3,8); 'DNF', 'savannas', changeTab(3,9); 'DNF', 'grassland', changeTab(3,10);
         'DBF', 'enf', changeTab(4,1); 'DBF', 'dnf', changeTab(4,3); 'DBF', 'mf', changeTab(4,5); 'DBF', 'woody savannas', changeTab(4,8); 'DBF', 'savannas', changeTab(4,9); 'DBF', 'grassland', changeTab(4,10); 'DBF', 'croplands', changeTab(4,12); 'DBF', 'cropland/natural vegetation mosaics', changeTab(4,14);
         'MF', 'enf', changeTab(5,1); 'MF', 'dnf', changeTab(5,3); 'MF', 'dbf', changeTab(5,4); 'MF', 'woody savannas', changeTab(5,8); 'MF', 'savannas', changeTab(5,9); 'MF', 'grassland', changeTab(5,10); 'MF', 'croplands', changeTab(5,12);
         'Woody Savannas', 'enf', changeTab(8,1); 'Woody Savannas', 'dnf', changeTab(5,3); 'Woody Savannas', 'dbf', changeTab(5,4); 'Woody Savannas', 'mf', changeTab(5,5); 'Woody Savannas', 'savannas', changeTab(5,9); 'Woody Savannas', 'grassland', changeTab(5,10); 'Woody Savannas', 'croplands', changeTab(5,12); 'Woody Savannas', 'cropland/natural vegetation mosaics', changeTab(5,14);
         'Savannas', 'enf', changeTab(9,1); 'Savannas', 'dbf', changeTab(9,4); 'Savannas', 'mf',changeTab(9,5); 'Savannas', 'woody savannas', changeTab(9,8); 'Savannas', 'grassland', changeTab(9,10); 'Savannas', 'croplands', changeTab(9,12); 'Savannas', 'cropland/natural vegetation mosaics', changeTab(9,14);
         'Grassland', 'enf', changeTab(10,1); 'Grassland', 'dnf', changeTab(10,3); 'Grassland', 'dbf', changeTab(10,4); 'Grassland', 'mf', changeTab(10,5); 'Grassland', 'woody savannas', changeTab(10,8); 'Grassland', 'savannas', changeTab(10,9); 'Grassland', 'croplands', changeTab(10,12); 'Grassland', 'cropland/natural vegetation mosaics', changeTab(10,14);
         'Croplands', 'dbf', changeTab(12,4); 'Croplands', 'mf', changeTab(12,5); 'Croplands', 'woody savannas', changeTab(12,8); 'Croplands', 'savannas', changeTab(12,9); 'Croplands', 'grassland', changeTab(12,10); 'Croplands', 'cropland/natural vegetation mosaics', changeTab(12,14);
         'Cropland/Natural Vegetation Mosaics', 'dbf', changeTab(14,4); 'Cropland/Natural Vegetation Mosaics', 'mf', changeTab(14,5); 'Cropland/Natural Vegetation Mosaics', 'woody savannas', changeTab(14,8); 'Cropland/Natural Vegetation Mosaics', 'savannas', changeTab(14,9); 'Cropland/Natural Vegetation Mosaics', 'grassland', changeTab(14,10); 'Cropland/Natural Vegetation Mosaics', 'croplands', changeTab(14,12)};

%% Figure plotting

f = figure('Position',[0 0 600 600]);


% figure a

ha = axes('Position',[-0.044 0.413 0.507 0.578]);

cmap = [57 210 252; 227 145 54]./255;

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

patch([-0.75 -0.65 -0.65 -0.75],[-0.95 -0.95 -0.89 -0.89],cmap(1,:))
patch([-0.75 -0.65 -0.65 -0.75],[-1.05 -1.05 -0.99 -0.99],cmap(2,:))

annotation(f,'textbox',[0.059 0.526 0.326 0.048],...
    'String',{'no land cover type change'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.059 0.498 0.291 0.048],...
    'String',{'land cover type change'},'LineStyle','none','FontSize',11);


% figure b

hb = axes('Position',[0.56 0.56 0.321 0.371]);

SK = SSankey(links(:,1), links(:,2), links(:,3));
SK.ColorList=[29 130 64; 
    157 192 66; 
    156 228 126; 
    125 195 125; 
    228 224 133; 
    223 225 0; 
    225 191 127; 
    226 226 130; 
    162 159 66;
    157 192 66; 
    156 228 126;
    125 195 125; 
    228 224 133; 
    223 225 0; 
    225 191 127; 
    29 130 64; 
    226 226 130;  
    162 159 66]./255;
SK.LabelLocation='left';
SK.draw()

annotation(f,'textbox',[0.526 0.928 0.099 0.051],...
    'String',{'2001'},'LineStyle','none','FontSize',12);
annotation(f,'textbox',[0.832 0.927 0.099 0.051],...
    'String',{'2020'},'LineStyle','none','FontSize',12);

% figure c

indl = find((K==1)&(FC_avg>0.05)&(FC_diff<0)&(MASKt==1)&(MASK_nochange==1));
indg = find((K==1)&(FC_avg>0.05)&(FC_diff>0)&(MASKt==1)&(MASK_nochange==1));

dataX1 = FC_diff(indl);
dataY1 = Tfcc(indl);
[xData1, yData1] = prepareCurveData( dataX1, dataY1 );
[fitresult1, gof1] = fit( xData1, yData1, 'poly1' );

dataX2 = FC_diff(indg);
dataY2 = Tfcc(indg);
[xData2, yData2] = prepareCurveData( dataX2, dataY2 );
[fitresult2, gof2] = fit( xData2, yData2, 'poly1' );

cmap = [182,110,151; 124,134,65]./255;

x=[1,2]; 
y = [fitresult1.p1,fitresult2.p1];

err3 = [0.0171, 0.0167];

hc = axes('Position',[0.147 0.070 0.388 0.375]);

b3 = bar(x(1),y(1),0.5);
hold on 
er3 = errorbar(x(1),y(1),err3(1));
b3.FaceColor = cmap(1,:);
b3.FaceAlpha = 0.7;
b3.EdgeColor = 'none';
b3.LineWidth = 1;
er3.LineWidth = 1.5;
er3.Color = cmap(1,:);
er3.CapSize = 12;

b4 = bar(x(2),y(2),0.5);
hold on 
er4 = errorbar(x(2),y(2),err3(2));
b4.FaceColor = cmap(2,:);
b4.FaceAlpha = 0.7;
b4.EdgeColor = 'none';
b4.LineWidth = 1.5;
er4.LineWidth = 1;
er4.Color = cmap(2,:);
er4.CapSize = 12;

ylim([-0.6,0])

ylabel({'Regional sensitivity ','of \DeltaTAC to \DeltaFC'});
set(hc,'FontSize',12,'LineWidth',1,'XTick',zeros(1,0),'YTickLabel',...
    {'-6','-5','-4','-3','-2','-1','0'},'YDir', 'reverse','Box','off');
title({'no land cover type change'},'FontSize',12);

annotation(f,'textbox',[0.196 0 0.158 0.081],...
    'String',{'Forest loss\newline    pixels'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.341 0 0.161 0.081],...
    'String',{'Forest gain\newline    pixels'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.143 0.407 0.123 0.059],...
    'String',{'(×10^-^3)'},'LineStyle','none','FontSize',12); % convertion into sensitivty of ΔTAC to 1% forest cover change (x10^-2)

%% 

% figure d

indl = find((K==1)&(FC_avg>0.05)&(FC_diff<0)&(MASKt==1)&(MASK_nochange==0));
indg = find((K==1)&(FC_avg>0.05)&(FC_diff>0)&(MASKt==1)&(MASK_nochange==0));

dataX1 = FC_diff(indl);
dataY1 = Tfcc(indl);
[xData1, yData1] = prepareCurveData( dataX1, dataY1 );
[fitresult1, gof1] = fit( xData1, yData1, 'poly1' );

dataX2 = FC_diff(indg);
dataY2 = Tfcc(indg);
[xData2, yData2] = prepareCurveData( dataX2, dataY2 );
[fitresult2, gof2] = fit( xData2, yData2, 'poly1' );

x=[1,2]; 
y = [fitresult1.p1,fitresult2.p1];

err3 = [0.0265, 0.0203];

hd = axes('Position',[0.568 0.070 0.388 0.375]);

b3 = bar(x(1),y(1),0.5);
hold on 
er3 = errorbar(x(1),y(1),err3(1));
b3.FaceColor = cmap(1,:);
b3.FaceAlpha = 0.7;
b3.EdgeColor = 'none';
b3.LineWidth = 1;
er3.LineWidth = 1.5;
er3.Color = cmap(1,:);
er3.CapSize = 12;

b4 = bar(x(2),y(2),0.5);
hold on 
er4 = errorbar(x(2),y(2),err3(2));
b4.FaceColor = cmap(2,:);
b4.FaceAlpha = 0.7;
b4.EdgeColor = 'none';
b4.LineWidth = 1.5;
er4.LineWidth = 1.5;
er4.Color = cmap(2,:);
er4.CapSize = 12;

ylim([-0.6,0])

title({'land cover type change'},'FontSize',12);
set(hd,'FontSize',12,'LineWidth',1,'XTick',zeros(1,0),'YTickLabel',...
    {'-6','-5','-4','-3','-2','-1','0'},'YDir', 'reverse','box','off');

annotation(f,'textbox',[0.623 0 0.158 0.081], ...
    'String',{'Forest loss\newline    pixels'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.768 0 0.161 0.081],...
    'String',{'Forest gain\newline    pixels'},'LineStyle','none','FontSize',11);
annotation(f,'textbox',[0.565 0.409 0.123 0.059],...
    'String',{'(×10^-^3)'},'LineStyle','none','FontSize',12); % convertion into sensitivty of ΔTAC to 1% forest cover change (x10^-2)

