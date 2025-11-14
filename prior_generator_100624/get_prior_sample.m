function [ms,ns,t_hist] = get_prior_sample(DISTS,z_vec,Nreals,do_stats)

if nargin < 4
    do_stats = 0;
end


% Function to get prior sample
Nz = length(z_vec);


% Prelocate vectors
ms = zeros(round(1.1*Nreals),Nz); % Lithologi vector
ns = zeros(round(1.1*Nreals),Nz); % Resistivity vector


% To report to user if issues arose during simulation
flag_vector = [0; 0];


% Time running time
tic

t_hist = cell(1,numel(DISTS.TYPES.types));

for i = 1:round(1.1*Nreals)

    % Get priors
    if do_stats == 1
        [m,layer_index,flag_vector,t_hist] = prior_m_reals_stats(DISTS,z_vec,flag_vector,t_hist);
    elseif do_stats == 0
        [m,layer_index,flag_vector] = prior_m_reals(DISTS,z_vec,flag_vector);
    end
    

    n = prior_n_reals(DISTS,m,layer_index);
    ms(i,:) = m;
    ns(i,:) = n;
    

    % Display progress
    if i/round(1.1*Nreals)*100 == round(i/round(1.1*Nreals)*100)
        clc
        disp('Generating priors')
        disp([num2str(i/round(1.1*Nreals)*100),' % Done'])
        disp(['Estimated time remaining: ', num2str(round(toc*(round(1.1*Nreals)/i-1))), ' seconds'])
    end    
end

fprintf('\nGenerated %.0f realizations\n', round(1.1*Nreals))
fprintf('Checking for identical realizations...\n')
[~, ims] = unique(ms, 'rows', 'stable');
fprintf('\nPruned %.0f dublicate realizations...\n', round(1.1*Nreals) - numel(ims))
ms = ms(ims, :);
ns = ns(ims, :);

while numel(ims) < Nreals
    fprintf('Generating %.0f more realizations and appending until all are unique...\n', round(1.1*Nreals - numel(ims)))
    ms_new = zeros(round(1.1*Nreals - numel(ims)), Nz);
    ns_new = zeros(round(1.1*Nreals - numel(ims)), Nz);

    for i = 1:round(1.1*Nreals - numel(ims))
        if do_stats == 1
            [m,layer_index,flag_vector,t_hist] = prior_m_reals_stats(DISTS,z_vec,flag_vector,t_hist);
        elseif do_stats == 0
            [m,layer_index,flag_vector] = prior_m_reals(DISTS,z_vec,flag_vector);
        end

        n = prior_n_reals(DISTS,m,layer_index);
        ms_new(i,:) = m;
        ns_new(i,:) = n;
    end

    ms = [ms; ms_new];
    ns = [ns; ns_new];
    
    [~, ims] = unique(ms, 'rows', 'stable');
    fprintf('\nPruned %.0f dublicate realizations...\n', round(1.1*Nreals - numel(ims)))
    ms = ms(ims, :);
    ns = ns(ims, :);

end

ms = ms(1:Nreals, :);
ns = ns(1:Nreals, :);

toc


% Provide warnings if issues occured
if flag_vector(1) == 1
    warning('Something went wrong and models might not represent your inputs. Consider if depths and thicknesses are reasonably chosen.')
end

if flag_vector(2) == 1
    warning('Somewhat succesfull. Number of layers possibly not a uniformly drawn')
end

if do_stats == 1
    figure; clf; set(gcf,'Color','w'); tiledlayout('flow');
    nexttile
    frequency_matrix = zeros(max(DISTS.TYPES.max_thick(1:end-1)),numel(DISTS.TYPES.types));
    edges = 0.5:1:max(DISTS.TYPES.max_thick(1:end-1))+0.5;
    for i = 1:numel(DISTS.TYPES.types)
        nexttile
        histogram(t_hist{i})
        ylabel('Thickness [m]')
        title(DISTS.TYPES.names(i))
        frequency_matrix(:,i) = histcounts(t_hist{i},edges)./Nreals;
    end
    figure; clf; set(gcf,'Color','w'); tiledlayout('horizontal');
    nexttile
    imagesc(frequency_matrix)
    set(gca,'xtick',DISTS.TYPES.types,'xticklabel',DISTS.TYPES.names)
    ylabel('Thickness [m]')
    colorbar
    set(gca,'Colormap',flipud(hot))
    nexttile
    n_layers = sum(diff(ms') ~= 0);
    edges = 0.5:1:max(n_layers)+0.5;
    histogram(n_layers,edges)
    xlabel('No. of layers')
    nexttile
    [~,~,~,counts] = count_category_all(ms');
    imagesc(counts)
    set(gca,'xtick',DISTS.TYPES.types,'xticklabel',DISTS.TYPES.names)
    ylabel('Depth [m]')
    nexttile
    % counts = additive_smoothing(counts,Nreals,0.01);
    E = zeros(Nz,1);
    for i = 1:Nz
        n_c = sum(counts(i,:) ~= 0);
        p = counts(i,counts(i,:) ~= 0)/Nreals;
        E(i) = -sum(p.*log(p)./log(5));
    end
    plot(E,z_vec,'-r')
    set(gca,'YDir','reverse')
    ylabel('Entropy')
    box off
end