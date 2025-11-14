close all; clear all; clc;
addpath data functions
rng(1)

Nreals = 1000000;
dmax = 46; 
dz = 1;

if ~exist(sprintf('prior_detailed_N%d_dmax%d.h5',Nreals,dmax))
    name = prior_generator_100624('prior_detailed',Nreals,dz,dmax);
else
    name = sprintf('prior_detailed_N%d_dmax%d.h5',Nreals,dmax);
end

fig_number = 97;

ns = h5read(name,'/M1')';
ms = h5read(name,'/M2')';
z_vec = h5readatt(name,'/M1','x')';
cmap = h5readatt(name,'/M2','cmap');
types = h5readatt(name,'/M2','class_name');
n_types = numel(types);

m_obs = [4*ones(1,8) 3*ones(1,7) 4*ones(1,2) 2*ones(1,3) 3*ones(1,3) 6*ones(1,22) 7*ones(1,1)];
d_obs = vector2matrix(m_obs,1:9,1);
split = (d_obs(24,3)+d_obs(24,6))/2;
d_obs(24,3) = split;
d_obs(24,6) = split;
n_layers_obs = 6;
t_layers_obs = [8 7 2 3 3.5 21.5 1];
t_types_obs = [4 3 4 2 3 6 7];


fig_number = summary_plots4(ms,name,fig_number,d_obs,n_layers_obs,t_layers_obs,t_types_obs,'Prior realizations');
set(gcf,'Position',[0 0 750 990])

nexttile(5,[4 4])
title('Prior marginal distribution, \rho(m)')

exportgraphics(gcf,'figures/fig3.png','Resolution',600)
