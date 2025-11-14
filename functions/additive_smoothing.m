function out = additive_smoothing(in,N,alpha)
out = (in+alpha)./(N+alpha*numel(in));
out = out./sum(out,2);
