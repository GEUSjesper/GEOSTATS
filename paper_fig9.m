close all; clear all; clc;
addpath functions data
rng(3)

Nreals = 1000000;
dmax = 60; 
dz = 1;

if ~exist(sprintf('prior_N%d_dmax%d.h5', Nreals, dmax))
    name = prior_generator_100624('prior', Nreals, dz, dmax);
else
    name = sprintf('prior_N%d_dmax%d.h5', Nreals, dmax);
end

ns = h5read(name, '/M1')';
ms = h5read(name, '/M2')';
z_vec = h5readatt(name, '/M1', 'x')';
cmap = h5readatt(name, '/M2', 'cmap');
types = h5readatt(name, '/M2', 'class_name');
n_types = numel(types);


%% Borehole locations and form

borehole_xs = [10, 20, 41];
borehole_depths = [50, 60, 50]./dz;
i_use = [4 5 6];

influence_dist = 20; % borehole influence cone exponent [m]

max_probability = 0.75;


%% Constructing borehole information
n_boreholes = numel(borehole_xs);
m_obs = ms(i_use, :);


% uniform distribution
uniform_probability = 1/n_types;
d_uni = zeros(numel(z_vec),n_types) + uniform_probability;


% profile x's
dx = 1;
x = dx:dx:50;
elev = 5 * sin(x / 4);
borehole_elevs = 5 * sin(borehole_xs / 4);


%% Sampling

% Prelocate various statistical matrices
modes_weights = nan(numel(z_vec), numel(x));
total_counts_weights = nan(1, numel(x));
max_counts_weights = nan(numel(z_vec), numel(x));
counts_all_weights = nan(numel(z_vec), numel(x), numel(types));
Ns_weights = zeros(1, numel(x));

modes_layers = nan(numel(z_vec), numel(x));
total_counts_layers = nan(1, numel(x));
max_counts_layers = nan(numel(z_vec), numel(x));
counts_all_layers = nan(numel(z_vec), numel(x), numel(types));
Ns_layers = zeros(1, numel(x));

% msConst = parallel.pool.Constant(ms);

tic

parfor i_x = 1:numel(x)

    % ms = msConst.Value;

    d_lik_weights = ones(1, Nreals); % prelocate likelihood vector
    d_lik_layers = ones(1, Nreals); % prelocate likelihood vector

    for i_b = 1:n_boreholes
        d_obs = d_uni;
        distance_to_borehole = abs(borehole_xs(i_b) - x(i_x));
        elevation_difference = round(elev(i_x) - borehole_elevs(i_b));
        if distance_to_borehole < influence_dist 
            highest_obs = interp1([0, influence_dist], [max_probability, uniform_probability], distance_to_borehole);

            m_obs_min = max([1, 1 - elevation_difference]);
            m_obs_max = min([dmax - elevation_difference, borehole_depths(i_b)]);

            d_obs_min = max([1, 1 + elevation_difference]);
            d_obs_max = min([borehole_depths(i_b) + elevation_difference, dmax]);

            d_obs_b = vector2matrix(m_obs(i_b, m_obs_min:m_obs_max), 1:n_types, highest_obs);

            d_obs(d_obs_min:d_obs_max, :) = d_obs(d_obs_min:d_obs_max, :) .* d_obs_b;
        end

        [w_vec,d_max] = weights(d_obs, z_vec, 0.2, 3);

        d_lik_weights = d_lik_weights .* d_likelihood_weights(ms, d_obs, d_max, w_vec)';
        d_lik_layers = d_lik_layers .* d_likelihood_layers(ms, d_obs);
        % for im = 1:Nreals
        %     d_lik_weights(im) = d_lik_weights(im) .* d_likelihood(ms(im,:), d_obs, d_max, w_vec);
        %     d_lik_layers(im) = d_lik_layers(im) .* d_likelihood_layer(ms(im, :), d_obs);
        % end
    end


    % find accepted realizations
    max_lik_weights = d_lik_weights./max(d_lik_weights);
    r = rand(size(max_lik_weights));
    accept_weights = find(max_lik_weights >= r);
    

    % Counting statistics
    if accept_weights > 0
        [modes_weights(:,i_x), max_counts_weights(:,i_x), total_counts_weights(i_x), counts_all_weights(:,i_x,:)] = count_category_all(ms(accept_weights,:)', 1:n_types);
    else
        warning('No matches')
    end
    Ns_weights(i_x) = numel(accept_weights);


    max_lik_layers = d_lik_layers./max(d_lik_layers);
    r = rand(size(max_lik_layers));
    accept_layers = find(max_lik_layers >= r);
    
    if accept_layers > 0
        [modes_layers(:,i_x), max_counts_layers(:,i_x), total_counts_layers(i_x), counts_all_layers(:,i_x,:)] = count_category_all(ms(accept_layers,:)', 1:n_types);
    else
        warning('No matches')
    end
    Ns_layers(i_x) = numel(accept_layers);

    if mod(i_x, 5) == 0 && labindex == 1
        fprintf('Progress: %.1f%%\n', 100*i_x/numel(x));
    end
end
toc


%%
fig_number = 97;
figure; clf; tiledlayout(4, 2, 'TileSpacing', 'compact')

ax1 = nexttile;
    hold on
    for i = 1:numel(x)
        img = imagesc(x(i), -dmax + elev(i) + 0.5 * dz, flipud(modes_weights(:,i)));
    end
    colormap(ax1, cmap);
    clim([0.5 n_types+0.5]);
    for i_b = 1:n_boreholes
        x_plot_image = [borehole_xs(i_b) - 0.3,  borehole_xs(i_b) + 0.3];
        y_plot_image = [borehole_elevs(i_b) - 0.5*dz, -borehole_depths(i_b) + borehole_elevs(i_b) + 0.5*dz];
        image('XData', x_plot_image, 'YData', y_plot_image, 'CData', repmat(m_obs(i_b, 1:borehole_depths(i_b))',1,1000));
        x_plot_line = [borehole_xs(i_b), borehole_xs(i_b), borehole_xs(i_b) + 0.6, borehole_xs(i_b) + 0.6, borehole_xs(i_b)] - 0.3;
        y_plot_line = [-borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, borehole_elevs(i_b) + 0.1, borehole_elevs(i_b) + 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1];
        plot(x_plot_line, y_plot_line,'-k')
    end
    plot(x,elev,'-k')
    xlim([0.5 max(x)+0.5])
    ylim([-65 10])
    title('Mode (Weighted approach)', 'FontWeight', 'normal')
    ylabel('Depth [m]')
    box on
    set(ax1, 'Layer', 'top');
    text(0.025,0.9,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','k'); fig_number = fig_number+1;

ax2 = nexttile;
    hold on
    for i = 1:numel(x)
        img = imagesc(x(i), -dmax + elev(i) + 0.5 * dz, flipud(modes_layers(:,i)));
    end
    colormap(ax2, cmap);
    clim([0.5 n_types+0.5]);
    col1 = colorbar;
    set(col1, 'YDir', 'reverse', 'TickLabels', types, 'Ticks', [linspace(1,n_types,n_types)]);
    for i_b = 1:n_boreholes
        x_plot_image = [borehole_xs(i_b) - 0.3,  borehole_xs(i_b) + 0.3];
        y_plot_image = [borehole_elevs(i_b) - 0.5*dz, -borehole_depths(i_b) + borehole_elevs(i_b) + 0.5*dz];
        image('XData', x_plot_image, 'YData', y_plot_image, 'CData', repmat(m_obs(i_b, 1:borehole_depths(i_b))',1,1000));
        x_plot_line = [borehole_xs(i_b), borehole_xs(i_b), borehole_xs(i_b) + 0.6, borehole_xs(i_b) + 0.6, borehole_xs(i_b)] - 0.3;
        y_plot_line = [-borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, borehole_elevs(i_b) + 0.1, borehole_elevs(i_b) + 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1];
        plot(x_plot_line, y_plot_line,'-k')
    end
    plot(x, elev, '-k')
    xlim([0.5 max(x)+0.5])
    ylim([-65 10])
    title('Mode (Grouping approach)', 'FontWeight', 'normal')
    ylabel('Depth [m]')
    box on
    set(ax2, 'Layer', 'top');
    text(0.025,0.9,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','k'); fig_number = fig_number+1;

ax3 = nexttile;
    hold on
    for i = 1:numel(x)
        E = log(1./(max_counts_weights(:,i)./total_counts_weights(i)));
        img = imagesc(x(i), -dmax + elev(i) + 0.5 * dz, flipud(E));
    end
    colormap(ax3, pink);
    clim([0 1]);
    for i_b = 1:n_boreholes
        x_left  = borehole_xs(i_b) - 0.3;
        x_right = borehole_xs(i_b) + 0.3;
        for k = 1:borehole_depths(i_b)
            y_top = borehole_elevs(i_b) - (k-1)*dz;
            y_bottom = borehole_elevs(i_b) - k*dz;
        
            x_patch = [x_left, x_right, x_right, x_left];
            y_patch = [y_top, y_top, y_bottom, y_bottom];
        
            patch(x_patch, y_patch, cmap( m_obs(i_b, k), :), 'EdgeColor', 'none');
        end
       
        x_plot_line = [borehole_xs(i_b), borehole_xs(i_b), borehole_xs(i_b) + 0.6, borehole_xs(i_b) + 0.6, borehole_xs(i_b)] - 0.3;
        y_plot_line = [-borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, borehole_elevs(i_b) + 0.1, borehole_elevs(i_b) + 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1];
        plot(x_plot_line, y_plot_line,'-k')
    end
    plot(x,elev,'-k')
    xlim([0.5 max(x)+0.5])
    ylim([-65 10])
    title('Entropy (Weighted approach)', 'FontWeight', 'normal')
    ylabel('Depth [m]')
    box on
    set(ax3, 'Layer', 'top');
    text(0.025,0.9,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','k'); fig_number = fig_number+1;

ax4 = nexttile;
    hold on
    for i = 1:numel(x)
        E = log(1./(max_counts_layers(:,i)./total_counts_layers(i)));
        img = imagesc(x(i), -dmax + elev(i) + 0.5 * dz, flipud(E));
    end
    colormap(ax4, pink);
    clim([0 1]);
    col1 = colorbar;
    ylabel(col1, 'Entropy')
    for i_b = 1:n_boreholes
        x_left  = borehole_xs(i_b) - 0.3;
        x_right = borehole_xs(i_b) + 0.3;
        for k = 1:borehole_depths(i_b)
            y_top = borehole_elevs(i_b) - (k-1)*dz;
            y_bottom = borehole_elevs(i_b) - k*dz;
        
            x_patch = [x_left, x_right, x_right, x_left];
            y_patch = [y_top, y_top, y_bottom, y_bottom];
        
            patch(x_patch, y_patch, cmap( m_obs(i_b, k), :), 'EdgeColor', 'none');
        end
        x_plot_line = [borehole_xs(i_b), borehole_xs(i_b), borehole_xs(i_b) + 0.6, borehole_xs(i_b) + 0.6, borehole_xs(i_b)] - 0.3;
        y_plot_line = [-borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, borehole_elevs(i_b) + 0.1, borehole_elevs(i_b) + 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1];
        plot(x_plot_line, y_plot_line,'-k')
    end
    plot(x,elev,'-k')
    xlim([0.5 max(x)+0.5])
    ylim([-65 10])
    title('Entropy (Grouping approach)', 'FontWeight', 'normal')
    ylabel('Depth [m]')
    box on
    set(ax4, 'Layer', 'top');
    text(0.025,0.9,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','k'); fig_number = fig_number+1;


lith = 2;
ax5 = nexttile;
    hold on
    for i = 1:numel(x)
        P = squeeze(counts_all_weights(:,i,lith))./total_counts_weights(i);
        img = imagesc(x(i), -dmax + elev(i) + 0.5 * dz, flipud(P));
    end
    colormap(ax5, P_color);
    clim([0 1]);
    for i_b = 1:n_boreholes
        x_left  = borehole_xs(i_b) - 0.3;
        x_right = borehole_xs(i_b) + 0.3;
        for k = 1:borehole_depths(i_b)
            y_top = borehole_elevs(i_b) - (k-1)*dz;
            y_bottom = borehole_elevs(i_b) - k*dz;
        
            x_patch = [x_left, x_right, x_right, x_left];
            y_patch = [y_top, y_top, y_bottom, y_bottom];
        
            patch(x_patch, y_patch, cmap( m_obs(i_b, k), :), 'EdgeColor', 'none');
        end
        x_plot_line = [borehole_xs(i_b), borehole_xs(i_b), borehole_xs(i_b) + 0.6, borehole_xs(i_b) + 0.6, borehole_xs(i_b)] - 0.3;
        y_plot_line = [-borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, borehole_elevs(i_b) + 0.1, borehole_elevs(i_b) + 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1];
        plot(x_plot_line, y_plot_line,'-k')
    end
    plot(x,elev,'-k')
    xlim([0.5 max(x)+0.5])
    ylim([-65 10])
    title('Probability of Till (Weighted approach)','FontWeight','normal')
    ylabel('Depth [m]')
    box on
    set(ax5, 'Layer', 'top');
    text(0.025,0.9,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','k'); fig_number = fig_number+1;

ax6 = nexttile;
    hold on
    for i = 1:numel(x)
        P = squeeze(counts_all_layers(:,i,lith))./total_counts_layers(i);
        img = imagesc(x(i), -dmax + elev(i) + 0.5 * dz, flipud(P));
    end
    colormap(ax6, P_color(100));
    clim([0 1]);
    col1 = colorbar;
    ylabel(col1, 'Probability of Till')
    for i_b = 1:n_boreholes
        x_left  = borehole_xs(i_b) - 0.3;
        x_right = borehole_xs(i_b) + 0.3;
        for k = 1:borehole_depths(i_b)
            y_top = borehole_elevs(i_b) - (k-1)*dz;
            y_bottom = borehole_elevs(i_b) - k*dz;
        
            x_patch = [x_left, x_right, x_right, x_left];
            y_patch = [y_top, y_top, y_bottom, y_bottom];
        
            patch(x_patch, y_patch, cmap( m_obs(i_b, k), :), 'EdgeColor', 'none');
        end
        x_plot_line = [borehole_xs(i_b), borehole_xs(i_b), borehole_xs(i_b) + 0.6, borehole_xs(i_b) + 0.6, borehole_xs(i_b)] - 0.3;
        y_plot_line = [-borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, borehole_elevs(i_b) + 0.1, borehole_elevs(i_b) + 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1, -borehole_depths(i_b) + borehole_elevs(i_b) - 0.1];
        plot(x_plot_line, y_plot_line,'-k')
    end
    plot(x,elev,'-k')
    xlim([0.5 max(x)+0.5])
    ylim([-65 10])
    title('Probability of Till (Grouping approach)','FontWeight','normal')
    ylabel('Depth [m]')
    box on
    set(ax6, 'Layer', 'top');
    text(0.025,0.9,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','k'); fig_number = fig_number+1;
   

ax7 = nexttile;
    x_long = 0:0.1:53;
    first = max([uniform_probability*ones(size(x_long)); max_probability - (max_probability - uniform_probability) * abs(borehole_xs(1)-x_long) / influence_dist]);
    second = max([uniform_probability*ones(size(x_long)); max_probability - (max_probability - uniform_probability) * abs(borehole_xs(2)-x_long) / influence_dist]);
    third = max([uniform_probability*ones(size(x_long)); max_probability - (max_probability - uniform_probability) * abs(borehole_xs(3)-x_long) / influence_dist]);

    plot(x_long,first,'-k')
    hold on
    plot(x_long,second,'--k')
    plot(x_long,third,'-.k')
    xlim([0.5 max(x)+0.5])
    xlabel('Distance along profile')
    ylabel('Obs probability')
    title('Applied probability of observation', 'FontWeight', 'normal')
    legend('1st well','2nd well','3rd well', 'Location', 'northeast') 
    text(0.025,0.9,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','k'); fig_number = fig_number+1;


ax8 = nexttile;
    plot(x, Ns_weights, '-k')
    hold on
    plot(x, Ns_layers, '--k')
    xlim([0.5 max(x)+0.5])
    xlabel('Distance along profile')
    ylabel('Count')
    title('Number of accepted realizations', 'FontWeight', 'normal')
    legend('Weighted method', ' Grouping method', 'Location', 'southeast')   
    set(ax8, 'YScale','log')
    text(0.025,0.9,[char(fig_number),')'],'Units','normalized','FontSize',12,'Color','k'); fig_number = fig_number+1;


exportgraphics(gcf,'figures/fig9.png','Resolution',600)


function d_lik = d_likelihood_weights(m, d_obs, d_max, wvec)

    [Nreals, n_z] = size(m);
    n_class = size(d_obs, 2);

    d_obs_norm = d_obs ./ d_max;  
    d_obs_w = d_obs_norm .^ (wvec(:));  

    idx = sub2ind([n_z, n_class], repmat((1:n_z), Nreals, 1), m);

    vals = d_obs_w(idx);

    d_lik = prod(vals, 2); 

end



function d_lik = d_likelihood_layers(m, d_obs)

    [n_z, n_class] = size(d_obs);
    Nrealz = size(m,1);
    d_lik = ones(1, Nrealz);
    m_obs = matrix2vector(d_obs);
    shift = diff(m_obs) ~= 0;
    shifts = find(shift == 1);
    shifts = [0 shifts n_z];
          
    for i_m = 1:Nrealz
        
        m_d = zeros(n_z, 1);
        for i = 1:n_class
           m_d(m(i_m,:) == i) = d_obs(m(i_m,:) == i,i);
        end
        
        for is = 1:numel(shifts)-1
            d_obs_segment = m_d(shifts(is)+1:shifts(is+1),:);
            lik_segment = mean(d_obs_segment);
            d_lik(i_m) = d_lik(i_m).*lik_segment;
        end
    end 

end