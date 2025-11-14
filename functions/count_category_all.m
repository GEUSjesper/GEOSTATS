function [mode,max_count,N_obs,counts] = count_category_all(obs,types)
if nargin < 2
    % types = unique(obs);
    types = 1:max(obs,[],'all');
end

n_types = numel(types);
N_obs = size(obs,2);
N_ds = size(obs,1);

counts = zeros(N_ds,n_types);
for i = 1:N_ds
    for j = 1:n_types
        counts(i,j) = sum(obs(i,:) == types(j));
    end
end

[max_count,mode_idx] = max(counts,[],2);
mode = types(mode_idx);
