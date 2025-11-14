function [d_lik] = d_likelihood_sample(m,d_obs,c)

if nargin<3
    c = 15;
end

[n_z,n_class] = size(d_obs);
[d_max,m_obs] = max(d_obs,[],2);

layers = cumsum([1; diff(m_obs) ~= 0]);
layer_no = max(layers);
sample_no = min([n_z c]);
sample_no = max([sample_no layer_no]);

indeces = 1:n_z;
idx = zeros(1,sample_no);

for i = 1:layer_no
    if sum(layers == i) > 1
        idx(i) = randsample(indeces(layers == i),1);
    else
        idx(i) = indeces(layers == i);
    end
end

indeces(idx(1:max(layers))) = [];

idx(idx == 0) = randsample(indeces,sample_no-layer_no);

m_d = zeros(n_z,1);
for i = 1:n_class
   m_d(m==i) = d_obs(m==i,i)./d_max(m==i);
end

d_lik = prod(m_d(idx));
