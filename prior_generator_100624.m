% close all; clear all; clc; rng(1)

% 1D prior generator for INTEGRATE project
% Just change filename
% filename = 'prior';



%%
function name = prior_generator_100624(filename,Nreals,dz,dmax,doPlot)
% Nreals = 10000;
if nargin < 3
    dz = 1;
end
if nargin < 4
    dmax = 90;
end
if nargin < 5
    doPlot = 1;
end

% Add subfolders with code
addpath prior_generator_100624 data


% Read excel file
warning('OFF','MATLAB:table:ModifiedAndSavedVarnames')
T = readtable(filename);


% Assigns lithologies
include = ~isnan(T.No_);
types = T.Class(include)';
min_thick = T.MinThickness(include)';
max_thick = T.MaxThickness(include)';
res_means = T.Resistivity(include)';
res_unc = T.ResistivityUncertainty(include)';
cmap = [T.Red(include),T.Green(include),T.Blue(include)]./255;


% Assign sections
N_sections = sum(~isnan(T.Unit));
for i = 1:N_sections
    section_types{i} = eval(T.Classes{i});
end


% TYPES
DISTS.TYPES.names = types;
DISTS.TYPES.types = 1:numel(types);
DISTS.TYPES.min_thick = min_thick;
DISTS.TYPES.max_thick = max_thick;
DISTS.TYPES.res = res_means;
DISTS.TYPES.res_unc = res_unc;


% TOP
DISTS.SECTIONS.N_sections = N_sections; % number of sections
DISTS.SECTIONS.types = section_types; % layer types possible in prequaternary
DISTS.SECTIONS.min_layers = T.MinNo_OfLayers; % min number of layers
DISTS.SECTIONS.max_layers = T.MaxNo_OfLayers; % max number of layers
DISTS.SECTIONS.min_sections = T.MinUnitThickness; % minimum thickness of layers [m]
DISTS.SECTIONS.max_sections = T.MaxUnitThickness; % maximum thickness of layers [m]
DISTS.SECTIONS.frequency = T.Frequency; % maximum thickness of layers [m]


% Get prior samples
z_vec = 1:dz:dmax;
[ms,ns] = get_prior_sample(DISTS,z_vec,Nreals,0);


% Write hdf5 file
name = sprintf('%s_N%d_dmax%d.h5',filename,Nreals,dmax);

cd data
try;delete(name);end
h5create(name,'/M1',size(ns'))
h5writeatt(name,'/','Creation date', date);
h5write(name,'/M1',ns')
h5writeatt(name,'/M1','is_discrete', 0);
h5writeatt(name,'/M1','name', 'Resistivity');
h5writeatt(name,'/M1','x', z_vec-dz/2); 
h5writeatt(name,'/M1','clim', [.1 2600]);
cmap_res = flj_log();
h5writeatt(name,'/M1','cmap', cmap_res);

h5create(name,'/M2',size(ms'))
h5write(name,'/M2',ms')
h5writeatt(name,'/M2','is_discrete', 1);
h5writeatt(name,'/M2','name', 'lithology');
h5writeatt(name,'/M2','class_name', types);
h5writeatt(name,'/M2','class_id', 1:length(types));
h5writeatt(name,'/M2','clim', [1 length(types)]);
h5writeatt(name,'/M2','cmap', cmap);

cd ../
%% Plot figures
if doPlot == 1
% Resistivity figure
figure; clf; set(gcf,'Color','w'); tl = tiledlayout('flow','TileSpacing','compact'); title(tl,'Resistivity distributions','FontSize',24)
for i = 1:numel(types)
    nexttile
    x=-1:0.01:4;
    y=normpdf(x,log10(res_means(i)),res_unc(i)*log10(res_means(i)));
    res_plot(10.^x,y)
    title(types(i))
    xlabel('Resistivity [\Omegam]')
    set(gca,'XTick',[0.1 1 10 100 1000],'XTickLabels',num2str([0.1 1 10 100 1000]'))
end
c = colorbar;
c.Layout.Tile = 'south';
set(c,'Ticks',log10([0.1 1 3.2 10 32 100 320 1000 3200]),'TickLabels',num2str([0.1 1 3.2 10 32 100 320 1000 3200]'))


% Lithologies of prior
figure; clf; set(gcf,'Color','w'); tl = tiledlayout('vertical'); title(tl,'Prior realizations','FontSize',24);
sp(1) = nexttile;
imagesc(1:100,z_vec,ms(1:100,:)')
hold on
ylabel('depth [m]');
title('Lithostratigraphy')
colormap(gca,cmap)
clim([1,length(types)])
col1 = colorbar;
col1.Ticks = [linspace(1.33,length(types)-0.34,length(types))];
col1.TickLabels = types;
set(col1, 'YDir', 'reverse' );
xlabel('Real #')


% Resistivities of prior
sp(2) = nexttile;
imagesc(1:100,z_vec,ns(1:100,:)')
hold on
ylabel('depth [m]');
title('Resistivity')
cmap_res = flj_log();
colormap(gca,cmap_res)
clim([0.1 2600])
col1 = colorbar();
set(col1,'XTick',[0.1 0.3 1 3.2 10 32 100 316 1000 2600]);
set(gca,'Colormap',cmap_res)
title(col1,'Resistivity [\Omegam]','color','k')
xlabel('Real #')
set(gca,'ColorScale','log')


linkprop(sp,{'XLim','YLim'});

end

