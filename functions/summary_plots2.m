function [fig_number,dist] = summary_plots2(ms,name,fig_number,d_obs,n_layers_obs,t_layers_obs,t_types_obs,tit)


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
n_types = numel(types);
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
tiledlayout(3,6,'TileSpacing','tight');
set(gcf,'Color','w')


% Plot observation
if exist('d_obs','var')
    nexttile(1,[1 3])
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
    title('d_{obs}','FontSize',14)
    set(gca,'xtick',[],'FontSize',12)
    xtickangle(90)
    ylabel('Depth [m]')
    clim([0 1])
    text(0.025,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','w'); fig_number = fig_number+1;
end


% Plot 100 models
nexttile(8,[1 5])
imagesc(ms(1:min([Nr 100]),:)')
clim([0.5 n_types+0.5]); xlim([0.5 100.5]);
set(gca,'fontsize', 12) 
set(gca,'Colormap',cmap); clim([0.5 n_types+0.5]);
c = colorbar; set(c,'YDir','reverse','xtick',1:n_types,'xticklabel',types,'FontSize',12)
% xlabel('Real #')
title([tit,', N = ',num2str(Nr)])
text(0.01,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','w'); fig_number = fig_number+1;


% Plot marginal distribution
nexttile(4,[1 3])
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
title('Marginal distribution')
% set(gca,'xtick',1:n_types,'xticklabel',types,'FontSize',12)
xtickangle(75)
ylabel('Depth [m]','FontSize',12)
clim([0 1])
if exist('m_obs')
    hold on
    plot(m_obs,1:Nm,'.r','MarkerSize',10)
end
text(0.025,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','w'); fig_number = fig_number+1;


% Plot entropy
nexttile(7,[1 1])
plot(E,z_vec,'-b')
set(gca,'YDir','reverse','FontSize',12)
title('Entropy')
xlabel('Entropy','FontSize',12)
xlim([0 1])
ylim([0 Nm])
ylabel('Depth [m]','FontSize',12)
box off
text(0.1,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','k'); fig_number = fig_number+1;


% Plot number of layers
nexttile(13,[1 3])
histogram('BinCounts',layer_counts,'BinEdges',edges);
title('Number of layers')
box off
set(gca,'FontSize',12)
xlabel('Number of layers','FontSize',12)
ylabel('Realizations','FontSize',12)
if exist('n_layers_obs')
    try
    hold on
    plot(n_layers_obs,layer_counts(n_layers_obs),'.r','MarkerSize',50)
    catch
    end
end
text(0.025,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','k'); fig_number = fig_number+1;


% Plot thickness of layers
nexttile(16,[1 3])
imagesc(thickness_counts./sum(thickness_counts,1))
title('Layer thicknesses')
set(gca,'Colormap',flipud(bone))
cb = colorbar;
% set(gca,'xtick',1:n_types,'xticklabel',types,'FontSize',12)
xtickangle(75)
ylabel('Thickness [m]','FontSize',12)
ylabel(cb,'Frequency')
if exist('t_layers_obs') && exist('t_types_obs')
    hold on
    plot(t_types_obs,t_layers_obs,'.r','MarkerSize',20)
end
text(0.025,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','w'); fig_number = fig_number+1;

fprintf('Mean entropy for distribution: %f\n',mean(E,'all'))


