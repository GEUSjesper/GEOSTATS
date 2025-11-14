function [d_lik] = d_likelihood_weight(m,d_obs,zvec_obs,dz)

[n_z,n_class] = size(d_obs);
d_max = max(d_obs,[],2);


shifts = find(diff(m) ~= 0);    % find shifts in lithology
shifts = [-1 shifts' numel(zvec_obs)-1]+dz/2;   % add top and bottom boundary


a = 1./(log(n_z));
b = 3;


dists = abs(zvec_obs-shifts');
w = a*(1-max(1./(1+dists.^b),[],1));
% w = w./sum(w);


m_d = zeros(n_z,1);
for i = 1:n_class
   m_d(m==i) = d_obs(m==i,i)./d_max(m==i);
end


d_lik = prod(m_d.^(w'));

end