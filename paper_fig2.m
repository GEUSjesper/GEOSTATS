close all; clear all; clc;
addpath data functions
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

fig_number = 97;

m_obs = [2*ones(1,8) 1*ones(1,7) 2*ones(1,2) 1*ones(1,6) 3*ones(1,22) 4*ones(1,1)];
d_obs = vector2matrix(m_obs,1:5,1);
split = (d_obs(24,1)+d_obs(24,3))/2;
d_obs(24,1) = split;
d_obs(24,3) = split;
n_layers_obs = 6;
t_layers_obs = [8 7 2 6.5 21.5 1];
t_types_obs = [2 1 2 1 3 4];

fig_number = summary_plots3(ms,name,fig_number,d_obs,n_layers_obs,t_layers_obs,t_types_obs,'Prior realizations');
set(gcf,'Position',[0 0 750 990])

nexttile(5,[4 4])
title('Prior marginal distribution, \rho(m)')

exportgraphics(gcf,'figures/fig2.png','Resolution',600)

