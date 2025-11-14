function [fig_number,dist] = summary_plots3(ms,name,fig_number,d_obs,n_layers_obs,t_layers_obs,t_types_obs,tit)


% Start subplots from a) if unspecified
if nargin < 3
    fig_number = 97;
end

if nargin < 8
    tit = 'Realizations';
end

% Load prior information and calculate other relevant values
z_vec = h5readatt(name,'/M1','x')';
cmap = h5readatt(name,'/M2','cmap');
types = h5readatt(name,'/M2','class_name');
clear types
row1 = {'Meltwater' 'Till' 'Miocene' 'Miocene' 'Paleogene'};
row2 = {'deposits' '' 'sand' 'clay' 'clay'};
labelArray = [row1; row2]; 
types = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
% n_types = numel(types);
n_types = 5;
[Nr,Nm] = size(ms);

if exist('d_obs','var')
    if size(d_obs,1) < length(z_vec)
        z_vec = z_vec(1:size(d_obs,1));
        Nm = size(d_obs,1);
    end
end


% Do calculations
[counts,mode,E,layer_counts,thickness_counts,edges] = distribution_stats(ms,name);


% Make new figure
figure; clf
tiledlayout(10,8,'TileSpacing','tight');
set(gcf,'Color','w')


% Plot observation
if exist('d_obs','var')
    nexttile(1,[4 4])
    imagesc(d_obs)
    hold on
    ax = gca;
    for i = 0.5:1:size(d_obs,2)
        plot([i i],ax.YLim,'-k','Color',[0.3 0.3 0.3])
    end
    for i = 0.5:1:size(d_obs,1)
        plot(ax.XLim,[i i],'-k','Color',[0.3 0.3 0.3])
    end
    colorbar
    set(gca,'xtick',1:n_types,'xticklabel',types)
    xtickangle(-90)
    ylabel('Depth [m]')
    set(gca,'fontsize', 10) 
    ax.XAxis.FontSize = 9;
    ax.YAxis.FontSize = 10;
    clim([0 1])
    title('Observation, d_{obs}','FontSize',11,'FontWeight','normal')
    text(0.025,0.95,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','w'); fig_number = fig_number+1;
end


% Plot marginal distribution
nexttile(5,[4 4])
dist = counts./Nr;
imagesc(dist)
hold on
ax = gca;
for i = 0.5:1:size(d_obs,2)
    plot([i i],ax.YLim,'-k','Color',[0.3 0.3 0.3])
end
for i = 0.5:1:size(d_obs,1)
    plot(ax.XLim,[i i],'-k','Color',[0.3 0.3 0.3])
end
colorbar
set(gca,'xtick',1:n_types,'xticklabel',types)
xtickangle(-90)
ylabel('Depth [m]')
clim([0 1])
set(gca,'fontsize', 10) 
ax.XAxis.FontSize = 9;
ax.YAxis.FontSize = 10;
title('Marginal distribution','FontSize',11,'FontWeight','normal')
text(0.025,0.95,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','w'); fig_number = fig_number+1;


% Plot entropy
nexttile(33,[3 2])
plot(E,z_vec,'-b')
set(gca,'YDir','reverse','YTick',10:10:45,'FontSize',10)
xlabel('Entropy','FontSize',10)
ylabel('Depth [m]','FontSize',10)
xlim([0 1])
ylim([0 Nm])
set(gca,'fontsize', 10) 
box on
title('Entropy','FontSize',11,'FontWeight','normal')
text(0.1,0.95,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','k'); fig_number = fig_number+1;


% Plot 100 models
nexttile(35,[3 6])
imagesc(ms(1:min([Nr 100]),:)')
clim([0.5 n_types+0.5]); xlim([0.5 100.5]);
set(gca,'Colormap',cmap); clim([0.5 n_types+0.5]);
c = colorbar;
set(c,'YDir','reverse','xtick',1:n_types,'xticklabel',types,'FontSize',10)
xlabel('Realization number')
ylabel('Depth [m]','FontSize',10)
set(gca,'fontsize', 10) 
title([tit,', N = ',num2str(Nr)],'FontSize',11,'FontWeight','normal')
text(0.01,0.95,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','w'); fig_number = fig_number+1;


% Plot number of layers
nexttile(57,[3 4])
histogram('BinCounts',layer_counts,'BinEdges',edges);
box on
xlabel('Number of layers','FontSize',10)
ylabel('Realizations','FontSize',10)
if exist('n_layers_obs')
    try
    hold on
    plot(n_layers_obs,layer_counts(n_layers_obs),'.r','MarkerSize',25)
    catch
    end
end
set(gca,'fontsize', 10) 
title('Layer counts','FontSize',11,'FontWeight','normal')
text(0.025,0.95,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','k'); fig_number = fig_number+1;


% Plot thickness of layers
nexttile(61,[3 4])
imagesc(thickness_counts./sum(thickness_counts,1))
set(gca,'Colormap',flipud(bone))
cb = colorbar;
set(gca,'xtick',1:n_types,'xticklabel',types,'FontSize',9)
xtickangle(-90)
ylabel('Thickness [m]','FontSize',10)
set(gca,'fontsize', 10) 
ax.XAxis.FontSize = 9;
ax.YAxis.FontSize = 10;
ylabel(cb,'Frequency')
set(cb,'FontSize',10)
if exist('t_layers_obs') && exist('t_types_obs')
    hold on
    plot(t_types_obs,t_layers_obs,'.r','MarkerSize',20)
end
title('Layer thicknesses','FontSize',11,'FontWeight','normal')
text(0.025,0.95,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','w'); fig_number = fig_number+1;

fprintf('Mean entropy for distribution: %f\n',mean(E,'all'))


