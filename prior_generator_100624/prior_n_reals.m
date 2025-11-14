function n = prior_n_reals(INFO,m,layer_index)

% Initialize n vector
n = m;

% Input resistivities for each layer (add with 5% uncertainty)
for i = 1:layer_index
    for j = 1:max(INFO.TYPES.types)
        n(m == j & layer_index == i) = 10.^(log10(INFO.TYPES.res(j)) + INFO.TYPES.res_unc(j)*log10(INFO.TYPES.res(j))*randn);
    end
end

