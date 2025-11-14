close all; clear all; clc;

addpath functions data
rng(1)

Nreals = 100000;
dmax = 46; 
dz = 1;

if ~exist(sprintf('prior_N%d_dmax%d.h5',Nreals,dmax))
    name = prior_generator_100624('prior',Nreals,dz,dmax);
else
    name = sprintf('prior_N%d_dmax%d.h5',Nreals,dmax);
end

% if ~exist(sprintf('prior_detailed_N%d_dmax%d.h5',Nreals,dmax))
%     name = prior_generator_100624('prior_detailed',Nreals,dz,dmax);
% else
%     name = sprintf('prior_detailed_N%d_dmax%d.h5',Nreals,dmax);
% end

ns = h5read(name,'/M1')';
ms = h5read(name,'/M2')';
z_vec = h5readatt(name,'/M1','x')';
cmap = h5readatt(name,'/M2','cmap');
types = h5readatt(name,'/M2','class_name');
n_types = numel(types);


n_infos = 1:5:46;
certainty = [0.99 0.75 0.5];
n_info = numel(n_infos);
n_certainty = numel(certainty);
n_times = 10;


N = zeros(n_times,n_info,n_certainty);
d_lik = zeros(1,Nreals);
d_lik_real = zeros(1,Nreals);
N_real = zeros(n_times,min([sum(n_infos <= 46) n_info]),n_certainty);

% m_obs_real = [2*ones(1,8) 1*ones(1,7) 2*ones(1,2) 1*ones(1,6) 3*ones(1,22) 4*ones(1,1)];
m_obs_real = [4*ones(1,8) 3*ones(1,7) 4*ones(1,2) 2*ones(1,3) 3*ones(1,3) 6*ones(1,22) 7*ones(1,1)];

tic
for i_n = 1:n_info

    for i_obs = 1:n_times
        m_obs = ms(i_obs,:);


        for i_c = 1:n_certainty


            % Constructing borehole information
            d_obs = vector2matrix(m_obs,1:n_types,certainty(i_c));
            
            d_obs_real_temp = vector2matrix(m_obs_real,1:n_types,certainty(i_c));
            % split = (d_obs_real_temp(24,1)+d_obs_real_temp(24,3))/2;
            % d_obs_real_temp(24,1) = split;
            % d_obs_real_temp(24,3) = split;
            split = (d_obs_real_temp(24,3)+d_obs_real_temp(24,6))/2;
            d_obs_real_temp(24,3) = split;
            d_obs_real_temp(24,6) = split;

            parfor im = 1:Nreals
                d_lik(im) = d_likelihood_test(ms(im,:),d_obs,n_infos(i_n)); % Likelihood
                d_lik_real(im) = d_likelihood_test(ms(im,:),d_obs_real_temp,n_infos(i_n)); % Likelihood
            end
            d_lik = d_lik./max(d_lik);
            d_lik_real = d_lik_real./max(d_lik_real);


            % find accepted realizations
            r = rand(1,Nreals);
            accept = find(d_lik >= r);
            N(i_obs,i_n,i_c) = numel(accept);
            accept_real = find(d_lik_real >= r);
            N_real(i_obs,i_n,i_c) = numel(accept_real);
        end
    end
    loopprogress('Calculating likelihoods',i_n,n_info) 
end
toc

%% PLOT SUMMARY
figure; clf; set(gcf,'Color','w');
tiledlayout('horizontal','TileSpacing','compact');
set(gcf,'Position',[680 394 599 484])


clrs = [0.9 0.2 0.1;
    0.3 0.7 0.3;
    0.1 0.2 0.9];


nexttile
hold on
% for i_n = 1:n_times
    for i_c = 1:n_certainty
        % % plot(n_infos,N(:,:,i_c),'-','Color',[clrs(i_c,:),0.1])
        h(i_c) = plot(n_infos,mean(N(:,:,i_c)),'-','Color',clrs(i_c,:),'LineWidth',2);
        plot(n_infos(n_infos <= 46),mean(N_real(:,:,i_c)),'--','Color',clrs(i_c,:));
    end
% end

xlabel('Number of model parameters')
ylabel('Number of posterior realizations')
set(gca,'YScale','log','FontSize',14)
grid on
% ylim([0 15000])
xlim([1 n_infos(n_info)])
set(findall(gcf,'type','text'),'FontSize',14)
legend(h,num2str(certainty'))


ax1 = gca;
y = ax1.YLim;
ytick = ax1.YTick;
yyaxis right
ax2 = gca;
set(ax2,'YScale','log')
set(ax2,'YLim',y)
ax2.YTickLabel = ytick./Nreals*100;
ylabel('Percentage of realizations (%)')

box on

% title('Simple prior (k = 5)')
% exportgraphics(gcf,'figures/fig3.png','Resolution',300)


