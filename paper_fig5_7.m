close all; clear all; clc;
addpath functions data
rng(1)


Nreals = 1000000;
dmax = 46; 
dz = 1;

if ~exist(sprintf('prior_N%d_dmax%d.h5',Nreals,dmax))
    name = prior_generator_100624('prior',Nreals,dz,dmax);
else
    name = sprintf('prior_N%d_dmax%d.h5',Nreals,dmax);
end

ns = h5read(name,'/M1')';
ms = h5read(name,'/M2')';
z_vec = h5readatt(name,'/M1','x')';
cmap = h5readatt(name,'/M2','cmap');
types = h5readatt(name,'/M2','class_name');
n_types = numel(types);


%% Constructing borehole information

m_obs = [2*ones(1,8) 1*ones(1,7) 2*ones(1,2) 1*ones(1,6) 3*ones(1,22) 4*ones(1,1)];
d_obs = vector2matrix(m_obs,1:5,0.75);
split = (d_obs(24,1)+d_obs(24,3))/2;
d_obs(24,1) = split;
d_obs(24,3) = split;


%% Sampling Nm

d_lik = zeros(1,Nreals);


tic
parfor im = 1:Nreals
    d_lik(im) = d_likelihood_test(ms(im,:),d_obs);
    loopprogress('Calculating likelihoods',im,Nreals)
end


maxL = d_lik./max(d_lik);


% find accepted realizations
r = rand(1,Nreals);
accept = find(maxL >= r);


%% Figure X

fig_number = 97;

n_layers_obs = 6;
t_layers_obs = [8 7 2 6.5 21.5 1];
t_types_obs = [2 1 2 1 3 4];

[fig_number,dist] = summary_plots3(ms(accept,:),name,fig_number,d_obs,n_layers_obs,t_layers_obs,t_types_obs,'Posterior realizations');
set(gcf,'Position',[0 0 750 990])

nexttile(5,[4 4])
title('Marginal distribution, \sigma(m)')

fprintf('KL divergence between marginal distribution and obs: %f\n',KLdivergence(d_obs,dist))
[~,~,~,counts] = count_category_all(ms',1:n_types);
fprintf('KL divergence between marginal distribution and prior: %f\n',KLdivergence(counts./Nreals,dist))
fprintf('Ratio: %f\n',KLdivergence(d_obs,dist)./KLdivergence(counts./Nreals,dist))


%% Sampling Ng
t_layers_obs = [0 8 7 2 7 21 1];
t_layers_obs = cumsum(t_layers_obs);
d_obs_collapsed = vector2matrix(t_types_obs,1:5,0.75);
ms_collapsed = zeros(Nreals,n_layers_obs);
for i = 1:Nreals
    for j = 1:n_layers_obs
        idx = t_layers_obs(j)+1:t_layers_obs(j+1);
        ms_collapsed(i,j) = mode(ms(i,idx));
    end
end
d_lik = zeros(1,Nreals);


tic
parfor im = 1:Nreals
    d_lik(im) = d_likelihood_test(ms_collapsed(im,:),d_obs_collapsed);
    loopprogress('Calculating likelihoods',im,Nreals)
end


maxL = d_lik./max(d_lik);


% find accepted realizations
r = rand(1,Nreals);
accept = find(maxL >= r);


%% Figure 7

fig_number = 97;

n_layers_obs = 6;
t_layers_obs = [8 7 2 6.5 21.5 1];
t_types_obs = [2 1 2 1 3 4];

[fig_number,dist] = summary_plots3(ms(accept,:),name,fig_number,d_obs,n_layers_obs,t_layers_obs,t_types_obs,'Posterior realizations');
set(gcf,'Position',[0 0 750 990])
nexttile(1)
hold off
m_obs = [2*ones(1,16) 1*ones(1,14) 2*ones(1,4) 1*ones(1,13) 3*ones(1,43) 4*ones(1,2)];
d_obs2 = vector2matrix(m_obs,1:5,0.75);
imagesc(d_obs2)
hold on
ax = gca;
for i = 0.5:1:size(d_obs,2)
    plot([i i],ax.YLim,'-k','Color',[0.2 0.2 0.2])
end
for i = cumsum(t_layers_obs)*2+0.5
    plot(ax.XLim,[i i],'-k','Color',[0.2 0.2 0.2])
end
colorbar
title('Observation, d_{obs}','FontSize',11,'FontWeight','normal')
labelArray = {'Meltwater','Till','Miocene','Miocene','Paleogene'; 'deposits','','sand','clay','clay'}; 
xticks = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
set(gca,'XTick',1:5,'XTickLabel',xticks)
xtickangle(-90)
ylabel('Depth [m]')
ax.XAxis.FontSize = 9;
ax.YAxis.FontSize = 10;
clim([0 1])
text(0.025,0.95,'a)','Units','normalized','FontSize',12,'Color','w');

nexttile(5,[4 4])
title('Posterior marginal distribution, \sigma(m)')

fprintf('KL divergence between marginal distribution and obs: %f\n',KLdivergence(d_obs,dist))
[~,~,~,counts] = count_category_all(ms',1:n_types);
fprintf('KL divergence between marginal distribution and prior: %f\n',KLdivergence(counts./Nreals,dist))
fprintf('Ratio: %f\n',KLdivergence(d_obs,dist)./KLdivergence(counts./Nreals,dist))

exportgraphics(gcf,'figures/fig7.png','Resolution',600)

%% Calculate weights for temperature

z_layers_obs = cumsum(t_layers_obs)';
w = 1.5*T_weights(z_layers_obs,z_vec,2);

d_obs_w = d_obs;

for i = 1:numel(w)
    d_obs_w(i,:) = d_obs(i,:).^(1./w(i))./sum(d_obs(i,:).^(1./w(i)));
end


%% Sampling Nm with weights

d_lik = zeros(1,Nreals);

tic
parfor im = 1:Nreals
    d_lik(im) = d_likelihood_test(ms(im,:),d_obs_w);
    loopprogress('Calculating likelihoods',im,Nreals)
end

maxL = d_lik./max(d_lik);


% find accepted realizations
r = rand(1,Nreals);
accept = find(maxL >= r);


%% Figure 5

fig_number = 97;

n_layers_obs = 6;
t_layers_obs = [8 7 2 6.5 21.5 1];
t_types_obs = [2 1 2 1 3 4];

[fig_number,dist] = summary_plots3(ms(accept,:),name,fig_number,d_obs_w,n_layers_obs,t_layers_obs,t_types_obs,'Posterior realizations');
set(gcf,'Position',[0 0 750 990])

nexttile(5,[4 4])
title('Posterior marginal distribution, \sigma(m)')

fprintf('KL divergence between marginal distribution and obs: %f\n',KLdivergence(d_obs,dist))
[~,~,~,counts] = count_category_all(ms',1:n_types);
fprintf('KL divergence between marginal distribution and prior: %f\n',KLdivergence(counts./Nreals,dist))
fprintf('Ratio: %f\n',KLdivergence(d_obs,dist)./KLdivergence(counts./Nreals,dist))

exportgraphics(gcf,'figures/fig5.png','Resolution',600)

