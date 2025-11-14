function res_plot(x,y)
% This function makes colored resistivity plots
%
% x is the resistivity values
% y is counts, probability or similar
%
% JN, 21/11-2023

xs = logspace(log10(0.1),log10(2600),443);
cs = flj_log;

if numel(xs) > numel(x)
    x_new = logspace(log10(min(x)),log10(max(x)),443);
    y_new = interp1(x,y,x_new);
    x = x_new;
    y = y_new;
end

cs_new = zeros(length(x),3);
cs_new(:,1) = interp1(xs,cs(:,1),x);
cs_new(:,2) = interp1(xs,cs(:,2),x);
cs_new(:,3) = interp1(xs,cs(:,3),x);
cs_new(x > max(xs),:) = repmat(cs(end,:),sum(x > max(xs)),1);
cs_new(x < min(xs),:) = repmat(cs(1,:),sum(x < min(xs)),1); 

patch([x fliplr(x)],[y zeros(size(y))],log10([x fliplr(x)]))

colormap(cs_new)

set(gca,'Xscale','log')
set(gca,'Layer', 'top')
