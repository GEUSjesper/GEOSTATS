function cmap_P = P_color(n)

if nargin < 1
    n = 255;
end

cmap_P = zeros(n,3);
numbers = linspace(0,1,n);
cmap_short = [1 1 1;...
            0.7 0.7 0.7;...
            1 0 0;...
            1 0.6 0;...
            1 1 0;...
            0 1 0];
% for i = 1:10
    cmap_P(:,1) = interp1(0:0.2:1,cmap_short(:,1),numbers);
    cmap_P(:,2) = interp1(0:0.2:1,cmap_short(:,2),numbers);
    cmap_P(:,3) = interp1(0:0.2:1,cmap_short(:,3),numbers);
% end
