function [d_lik] = d_likelihood_test(m,d_obs,I)

% Likelihood function so that L = prod(d_obs(m(I)))
% m is realization vector of size n_z
% d_obs is probabilities of observation of size n_z,n_class
% I is indecies to use (default: all)

[n_z,n_class] = size(d_obs);
[d_max,m_obs] = max(d_obs,[],2);


if nargin<3
    I = n_z;
end

number_of_samples = min([n_z I]);

idx = randperm(n_z,number_of_samples);

m_d = zeros(n_z,1);
for i = 1:n_class
   m_d(m==i) = d_obs(m==i,i)./d_max(m==i);
end

d_lik = prod(m_d(idx));
