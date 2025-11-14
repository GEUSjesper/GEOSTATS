function [res_quant,N_obs] = count_quantile(obs,quantile_vector)

N_obs = size(obs,2);
N_ds = size(obs,1);

res_quant = zeros(N_ds,numel(quantile_vector));

for i = 1:N_ds
    res_quant(i,:) = quantile(obs(i,:),quantile_vector);
end
