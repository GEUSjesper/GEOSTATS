function [fig_number,dist] = summary_plots_pres(ms,name,fig_number,m_obs,n_layers_obs,t_layers_obs,t_types_obs,tit)


% Start subplots from a) if unspecified
if nargin < 3
    fig_number = 97;
end

if nargin < 8
    tit = 'Prior realizations';
end

% Load prior information and calculate other relevant values
z_vec = h5readatt(name,'/M1','x')';
cmap = h5readatt(name,'/M2','cmap');
types = h5readatt(name,'/M2','class_name');
n_types = numel(types);
[Nr,Nm] = size(ms);

if exist('m_obs')
    if length(m_obs) < length(z_vec)
        z_vec = z_vec(1:length(m_obs));
        Nm = length(m_obs);
    end
end

% Do calculations
[counts,mode,E,layer_counts,thickness_counts,edges] = distribution_stats(ms,name);


% Make new figure
figure; clf
tiledlayout(1,14,'TileSpacing','compact');
set(gcf,'Color','w')


% Plot observation
if exist('m_obs')
    nexttile(1)
    imagesc(m_obs')
    hold on
    if m_obs(24) == 6
        fill([0.5 1 1 0.5 0.5],[24.5 24.5 23.5 23.5 24.5],cmap(3,:),'LineStyle','none')
    else
        fill([0.5 1 1 0.5 0.5],[24.5 24.5 23.5 23.5 24.5],cmap(1,:),'LineStyle','none')
    end
    clim([0.5 n_types+0.5]); xlim([0.5 1.5]); ylabel('Depth [m]');
    title('Obs')
    set(gca,'xticklabels',[])
    set(gca,'fontsize', 12) 
    set(gca,'Colormap',cmap); clim([0.5 n_types+0.5]);
    text(0.05,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','w'); fig_number = fig_number+1;
else
    nexttile(1)
    plot(1,1)
    ylabel('Depth [m]');
    set(gca,'fontsize', 12) 
end


% Plot mode
nexttile(2)
imagesc(mode')
clim([0.5 n_types+0.5]); xlim([0.5 1.5]); ylabel('Depth [m]');
title('Mode')
set(gca,'xticklabels',[])
set(gca,'fontsize', 12) 
set(gca,'Colormap',cmap); clim([0.5 n_types+0.5]);
text(0.05,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','w'); fig_number = fig_number+1;


% Plot 100 models
nexttile(3,[1 5])
imagesc(ms(1:min([Nr 100]),:)')
clim([0.5 n_types+0.5]); xlim([0.5 100.5]);
set(gca,'fontsize', 12) 
set(gca,'Colormap',cmap); clim([0.5 n_types+0.5]);
% c = colorbar('southoutside'); set(c,'xtick',1:n_types,'xticklabel',types,'FontSize',12)
xlabel('Real #')
title([tit,', N = ',num2str(Nr)])
text(0.01,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','w'); fig_number = fig_number+1;


% Plot marginal distribution
nexttile(8,[1,2])
dist = counts./Nr;
imagesc(dist)
colorbar('southoutside')
title('Marginal distribution')
set(gca,'xtick',1:n_types,'xticklabel',types,'FontSize',12)
xtickangle(90)
ylabel('Depth [m]','FontSize',12)
clim([0 1])
if exist('m_obs')
    hold on
    plot(m_obs,1:Nm,'.r','MarkerSize',10)
end
text(0.025,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','k'); fig_number = fig_number+1;


% Plot entropy
nexttile(10)
plot(E,z_vec,'-b')
set(gca,'YDir','reverse','FontSize',12)
title('Entropy')
xlabel('Entropy','FontSize',12)
xlim([0 1])
ylim([0 Nm])
ylabel('Depth [m]','FontSize',12)
box off
text(0.05,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','k'); fig_number = fig_number+1;


% Plot number of layers
nexttile(11,[1,2])
histogram('BinCounts',layer_counts,'BinEdges',edges);
title('Number of layers')
box off
set(gca,'FontSize',12)
xlabel('Number of layers','FontSize',12)
ylabel('Realizations','FontSize',12)
if exist('n_layers_obs')
    hold on
    plot(n_layers_obs,layer_counts(n_layers_obs),'.r','MarkerSize',50)
end
text(0.025,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','k'); fig_number = fig_number+1;


% Plot thickness of layers
nexttile(13,[1,2])
imagesc(thickness_counts./sum(thickness_counts,1))
title('Layer thicknesses')
set(gca,'Colormap',flipud(bone))
colorbar('southoutside')
% ylabel(cb,'Frequency')
set(gca,'xtick',1:n_types,'xticklabel',types,'FontSize',12)
xtickangle(90)
ylabel('Thickness [m]','FontSize',12)
if exist('t_layers_obs') && exist('t_types_obs')
    hold on
    plot(t_types_obs,t_layers_obs,'.r','MarkerSize',20)
end
text(0.025,0.95,[char(fig_number),')'],'Units','normalized','FontSize',15,'Color','w'); fig_number = fig_number+1;

fprintf('Mean entropy for distribution: %f\n',mean(E,'all'))


