function D_OBS = boreholespread(d_obs, x_obs, y_obs, elev_obs, UTMX, UTMY, ELEVATION, influence_max, Nm, dz)


if nargin < 10
    dz = 1;
end

if nargin < 9 || isempty(Nm)
    Nm = size(d_obs,2);
end


% Number of boreholes and number of tTEM soundings
n_obs = numel(x_obs);
N = numel(UTMX);
n_types = size(d_obs,3);


% Calculate distances to the borehole
DISTS = zeros(n_obs,N);
for i = 1:n_obs
    DISTS(i,:) = sqrt((UTMX-x_obs(i)).^2+(UTMY-y_obs(i)).^2);
end


% An uniform observation
d_uni = ones(Nm,n_types);


% Prelocate observation matrix
D_OBS = 1 / n_obs * ones(N, Nm, n_types);
% D_OBS = 1/n_obs*zeros(N,Nm,n_types);
% D_OBS = zeros(N,Nm,3);


% Find depth of boreholes
borehole_depth = size(d_obs, 2) * dz;


% Loop N sounding locations
for i=1:N
    % Loop boreholes
    for i_obs = 1:n_obs
        % If within borehole influence, else use a uniform distribution
        if DISTS(i_obs,i) < influence_max
            d_obs_i = d_uni;
            d_obs_temp = squeeze(d_obs(i_obs,:,:));


            % z vector adjusted for elevation difference from observation and
            % borehole obersvation length
            elev_diff = round((ELEVATION(i)-elev_obs(i_obs))./dz);

    
            % Align boreholeinformation           
            if elev_diff > 0
                % d_obs_i(elev_diff:min([elev_diff+borehole_depth-1 end]),:) = d_obs_temp;
                d_obs_i(elev_diff:min([elev_diff+borehole_depth-1 Nm]),:) = d_obs_temp(1:min([borehole_depth Nm-elev_diff+1]),:);
            elseif elev_diff <= 0
                d_obs_i(1:borehole_depth + elev_diff, :) = d_obs_temp(abs(elev_diff) + 1:borehole_depth, :);
            end

    
            % Write observation with appropriate distance weighting
            % D_OBS(i,:,:) = squeeze(D_OBS(i,:,:)).*((d_obs_i*(1-DISTS(i_obs,i)./influence_max) + d_uni*DISTS(i_obs,i)./influence_max));
            % D_OBS(i,:,:) = squeeze(D_OBS(i,:,:)).*((d_obs_i*(1-(DISTS(i_obs,i)./influence_max).^2) + d_uni*(DISTS(i_obs,i)./influence_max).^2));
            % D_OBS(i,:,:) = squeeze(D_OBS(i,:,:))+(d_obs_i*(1-DISTS(i_obs,i)./influence_max) + d_uni*DISTS(i_obs,i)./influence_max);
            D_OBS(i,:,:) = squeeze(D_OBS(i,:,:))*(DISTS(i_obs,i)./influence_max).^2+d_obs_i*(1-DISTS(i_obs,i)./influence_max).^2;
        end
    end
    D_OBS(i,:,:) = D_OBS(i,:,:)./sum(D_OBS(i,:,:),3);
end

