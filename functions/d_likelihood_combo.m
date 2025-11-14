function [d_lik] = d_likelihood_combo(m,d_obs,n)

if nargin < 3
    n = 15;
end

[n_z,n_class] = size(d_obs);
[d_max,m_obs] = max(d_obs,[],2);

idx = randperm(n_z,n);

mat = zeros(n_z,n);
for i = 1:n
    mat(:,i) = abs((1:n_z)-idx(i));
end

[~,a] = min(mat,[],2);

m_d = zeros(n_z,1);
for i = 1:n_class
   m_d(m==i) = d_obs(m==i,i)./d_max(m==i);
end

lik = zeros(n,1);

for i = 1:n
    lik(i) = mean(m_d(a == i));
end

d_lik = prod(lik);


