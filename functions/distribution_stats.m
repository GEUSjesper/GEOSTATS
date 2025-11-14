function [counts,mode,E,layer_counts,thickness_counts,edges] = distribution_stats(ms,name)


types = h5readatt(name,'/M2','class_name');
n_types = numel(types);
[Nreals,Nz] = size(ms);


% Counts all classes
[mode,~,~,counts] = count_category_all(ms',1:n_types);


% Calculate entropy
E = zeros(Nz,1);
for i = 1:Nz
    p = counts(i,counts(i,:) ~= 0)/Nreals;
    E(i) = -sum(p.*log(p)./log(n_types));
end


% Calculate number of layers
n_layers = sum(diff(ms') ~= 0) + 1;
edges = 0.5:1:max(n_layers)+1.5;
layer_counts = histcounts(n_layers,edges);


% Calculate thickness matrix
thickness_counts = zeros(Nz,n_types);
ms_temp = [ms'; zeros(1,Nreals)]; % create ms with a 0 devider
ms_array = reshape(ms_temp,Nreals.*Nz+Nreals,1);
for i = 1:n_types
    counter = 0;
    for j = 1:numel(ms_array)
        if ms_array(j) == i
            counter = counter + 1;
        elseif ms_array(j) ~= i && counter > 0
            thickness_counts(counter,i) = thickness_counts(counter,i) + 1;
            counter = 0;
        end
    end
end


