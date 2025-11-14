function [d_lik] = d_likelihood_layer(m,d_obs)

m_obs = matrix2vector(d_obs);
shift = diff(m_obs) ~= 0;
shifts = find(shift == 1);
shifts = [0 shifts numel(m)];
[n_z,n_class] = size(d_obs);
d_max = max(d_obs,[],2);
d_lik_lay = 1;

% for is = 1:numel(shifts)-1
%     l_segment = numel(shifts(is)+1:shifts(is+1));
%     d_obs_segment = d_obs(shifts(is)+1:shifts(is+1),:);
%     m_segment = m(shifts(is)+1:shifts(is+1));
%     d_segment = vector2matrix(m_segment,1:n_class,1);
%     max_segment = max(d_obs_segment,[],2);
%     lik_segment = mean(d_obs_segment(logical(d_segment))./max_segment);
%     % lik_segment = prod(d_obs_segment(logical(d_segment))./max_segment).^(1/l_segment);
%     d_lik_lay = d_lik_lay.*lik_segment;
% end

m_d = zeros(n_z,1);
for i = 1:n_class
   m_d(m==i) = d_obs(m==i,i);%./d_max(m==i);
end

for is = 1:numel(shifts)-1
    d_obs_segment = m_d(shifts(is)+1:shifts(is+1),:);
    lik_segment = mean(d_obs_segment);
    d_lik_lay = d_lik_lay.*lik_segment;
end

% d_lik = d_lik_lay.^(1/sum(shifts));
d_lik = d_lik_lay;