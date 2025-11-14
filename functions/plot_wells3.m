function plot_wells3(ELEV,well_lps,well_liths,cmap_lith,well_width)
% To plot wells on top of plot_profile function

for i = 1:numel(well_liths)
    liths = well_liths{i};
    diffs = true;
    for j = 1:numel(liths)-1
        diffs = [diffs; ~(liths(j) == liths(j+1))];
    end
    liths = liths(diffs);
    depths = [find(diffs == 1)-1; numel(well_liths{i})];
    n_layers = numel(liths);
    for j = 1:n_layers
        x = [well_lps(i)-well_width well_lps(i)-well_width well_lps(i)+well_width well_lps(i)+well_width];
        y = ELEV(i) + [depths(j) depths(j+1) depths(j+1) depths(j)];
        patch(x,y,cmap_lith(liths(j),:));
    end
end
