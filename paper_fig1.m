close all; clear all; clc;
addpath data
rng(1)

cds= [255 85 0]/255;
cml = [219 182 0]/255;
cgl = [178 255 255]/255;
cgs = [231 251 251]/255;
types1 = {'Sand','Clay'};
types2 = {'Meltwater deposits','Till','Miocene sand','Miocene clay','Paleogene clay'};
types3 = {'Meltwater clay','Meltwater sand','Meltwater gravel',...
    'Clayley till','Sandy till','Miocene sand','Miocene clay','Paleogene clay','Organic rich deposit'};

depths = [0 8 15 17 20 23.5 45 46];
cs = [cml; cds; cml; cds; cds; cgs; cgl];

fig_number = 97;

figure; clf; set(gcf,'Color','w');
fig1 = tiledlayout(1,3,'TileSpacing','compact');
set(gcf,'Color','w','Position',[2092 109 1092 865]);

nexttile(1,[1 2])
for i = 1:numel(depths)-1
    patch([0 1 1 0 0],[depths(i) depths(i) depths(i+1) depths(i+1) depths(i)],cs(i,:),'LineWidth',1);
end
xlim([0 10])
ylim([0 46])
ylabel('Depth [m]')
set(gca,'YDir','reverse','XTick',[],'xcolor',[1 1 1],'FontSize',13)
text(1.05,0,'0','FontSize',13)
text(1.05,8,'8','FontSize',13)
text(1.05,15,'15','FontSize',13)
text(1.05,17,'17','FontSize',13)
text(1.05,20,'20','FontSize',13)
text(1.05,23.5,'23.5','FontSize',13)
text(1.05,45,'45','FontSize',13)
text(0.01,0.975,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','k'); fig_number = fig_number+1;


text(1.43,1,{'CLAY: sandy, slightly gravelly, yellow-brown, calcareous, "clayey till".'},'VerticalAlignment','baseline','FontSize',12)
text(1.43,9,{'SAND AND GRAVEL: unsorted, few stones, few clumps of clay,','yellow-brown, calcareous, "meltwater sand and gravel".'},'VerticalAlignment','baseline','FontSize',12)
text(1.43,15,{'CLAY: sandy, many slivers of sand, dark grey-brown, calcareous,','"clayey till".'},'VerticalAlignment','baseline','FontSize',12)
text(1.43,18.25,{'SAND AND GRAVEL: unsorted, yellow-brown, calcareous, "meltwater sand','and gravel". Sample taken from 20 m.'},'VerticalAlignment','baseline','FontSize',12)
text(1.43,21.5,{'GRAVEL AND STONES: very sandy, yellow-brown, calcareous, "meltwater','stones". Samples taken from 23 m.'},'VerticalAlignment','baseline','FontSize',12)
text(1.43,26,{'SAND: mostly fine, many slivers of clay, dark grey-brown, micareous,','glauconite, non-calcareous, "mica sand".'},'VerticalAlignment','baseline','FontSize',12)
text(1.43,45,{'CLAY: silty, slivers of sand, black, micareous, glauconite, non-calcareous,','"mica clay".'},'VerticalAlignment','baseline','FontSize',12)


nexttile
m_obs = [2*ones(1,8) 1*ones(1,7) 2*ones(1,2) 1*ones(1,6) 3*ones(1,22) 4*ones(1,1)];
d_obs = vector2matrix(m_obs,1:5,1);
split = (d_obs(24,1)+d_obs(24,3))/2;
d_obs(24,1) = split;
d_obs(24,3) = split;
imagesc(d_obs)
set(gca,'xtick',1:5,'xticklabel',types2,'FontSize',14)
xtickangle(30)
set(gca,'ytick',0.5:5:45.5,'yticklabel',num2str([0 5 10 15 20 25 30 35 40 45]'))
hold on
ax = gca;
for i = 0.5:1:size(d_obs,2)
    plot([i i],ax.YLim,'-k','Color',[0.3 0.3 0.3]./3)
end
for i = 0.5:1:size(d_obs,1)
    plot(ax.XLim,[i i],'-k','Color',[0.3 0.3 0.3]./3)
end


clim([0 1])
text(0.05,0.975,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','w'); fig_number = fig_number+1;
cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb,'Probability','FontSize',14)

a1 = gca;
xt = get(a1,'XTick');
a2 = copyobj(a1,fig1);
set(a2,'Color','none')
set(a2,'Ytick',[])
set(a2,'XAxisLocation','top')
set(a2,'XTickLabel',['c_1'; 'c_2'; 'c_3'; 'c_4'; 'c_5';])
xtickangle(a2,0)

exportgraphics(fig1,'figures/fig1.png','Resolution',600)

