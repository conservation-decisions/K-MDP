%Add dependencies
addpath('kmeans')
addpath('utils')
addpath('algorithms')
addpath('MDP_graph_representation')

%precision
p = 0.00000000000001;
%Discount factor
discount = 0.96;

%Load transition and reward matrices
load('problems/SONA/SONA_P_500_sim.mat');
load('problems/SONA/SONA_R_500_sim.mat');


P = Tr;
NS = size(P,1);

%Abstract state space size
K = [round(NS), round(NS/1.25), round(NS/1.5), round(NS/2), round(NS/4)];


tic;

%Perform value iteration
[Pol]=mdp_value_iteration(P, R, discount);%Get a policy using value iteration
%Evaluate the policy
[V,Q]= mdp_eval_policy_iterative_q(P, R, discount, Pol); % check that this is doing what it is supposed to; Get the value and Q values

time = toc;

Q(isnan(Q)) = 0;



%Initialize structures for statistics extraction

err_transitive = zeros(length(K),1);
t_transitive = zeros(length(K),1);
t_KMDP_transitive = zeros(length(K),1);

err_astar = zeros(length(K),1);
t_astar = zeros(length(K),1);
t_KMDP_astar = zeros(length(K),1);

err_kmeans = zeros(length(K),1);
t_kmeans = zeros(length(K),1);
t_KMDP_kmeans = zeros(length(K),1);


S2K_Qd = cell(length(K),1);
K2S_Qd = cell(length(K),1);

S2K_astar = cell(length(K),1);
K2S_astar = cell(length(K),1);


%RUN the algorithms


for i = 1:length(K)
    
  k = K(i);
   
 fprintf('Q_d\n');
 [PK_Qd, RK_Qd, S2K_Qd{i}, K2S_Qd{i}, PolicyK_Qd, PolKs_Qd, err_transitive(i), t_transitive(i), t_KMDP_transitive(i)] = QdKMDP(k, p, P, R, discount, Q, V);
                
 fprintf('a*_d\n');
 [PK_astar, RK_astar, S2K_astar, K2S_astar, PolicyK_astar, PolKs_astar, err_astar(i), t_astar(i), t_KMDP_astar(i)] = aStarKMDP(k, p, P, R, discount, V, Pol); 
    
 fprintf('K-means \n');
 [PK_kmeans, RK_Kmeans, S2K_kmeans, K2S_kmeans, PolicyK_kmeans, PolKs_kmeans, err_kmeans(i), t_kmeans(i), t_KMDP_kmeans(i)] = kmeansKMDP(k, P, R, discount, Q, V);

    
end


%****************************************************************************************************************************

%PLOT RESULTS: error, times and state space representations


%DATA to plot and extract

disp('Q_d')
err_transitive
t_transitive
t_KMDP_transitive

disp('astar')
err_astar
t_astar
t_KMDP_astar


disp('kmeans++')
err_kmeans
t_kmeans
t_KMDP_kmeans

%Store relevant data. Saves the last K-MDP

save('problems/SONA/results/S2K_kmeans.mat', 'S2K_kmeans');
save('problems/SONA/results/K2S_kmeans.mat', 'K2S_kmeans');
save('problems/SONA/results/S2K_Qd.mat', 'S2K_Qd');
save('problems/SONA/results/K2S_Qd.mat', 'K2S_Qd');
save('problems/SONA/results/S2K_astar.mat', 'S2K_astar');
save('problems/SONA/results/K2S_astar.mat', 'K2S_astar');


err_transitive = err_transitive * 100;
err_astar = err_astar * 100;
err_kmeans = err_kmeans * 100;

%Plot releveant data


            %Error

            figure;
            plot(K, err_transitive, 'r-o', 'LineWidth', 1)
            hold on;
            plot(K, err_astar, 'b-x', 'LineWidth',1)
            hold on;
            plot(K, err_kmeans, 'g-*', 'LineWidth',1)
            hold off;
            xlabel('K')
            ylabel('gap(%)');
            title('Sea Otter and Northern Abalone error');
            legend('Q^*_d K-MDP', 'Q^*_a K-MDP', 'k-means++ K-MDP');
            plot_name_fig = strcat('problems/SONA/results/gap.fig');
            plot_name_png = strcat('problems/SONA/results/gap.png');
            saveas(gcf, plot_name_fig);
            saveas(gcf, plot_name_png);
            
            
            %Time compute KDMP
            
            figure;
            plot(K, t_transitive, 'r-o', 'LineWidth', 1)
            hold on;
            plot(K, t_astar, 'b-x', 'LineWidth',1)
            hold on;
            plot(K, t_kmeans, 'g-*', 'LineWidth',1)
            hold off;
            xlabel('K')
            ylabel('time(sec.)');
            title('Sea Otter and Northern Abalone time to compute the K-MDP');
            legend('Q^*_d K-MDP', 'Q^*_a K-MDP', 'k-means++ K-MDP');
            plot_name_fig = strcat('problems/SONA/results/time-compute.fig');
            plot_name_png = strcat('problems/SONA/results/time-compute.png');
            saveas(gcf, plot_name_fig);
            saveas(gcf, plot_name_png);
            
            
            %Time to solve KDMP
            
            figure;
            plot(K, t_KMDP_transitive, 'r-o', 'LineWidth', 1)
            hold on;
            plot(K, t_KMDP_astar, 'b-x', 'LineWidth',1)
            hold on;
            plot(K, t_KMDP_kmeans, 'b-x', 'LineWidth',1)
            hold off;
            xlabel('K')
            ylabel('time(sec.)');
            title('Sea Otter and Northern Abalone time to solve the K-MDP');
            legend('Q^*_d K-MDP', 'Q^*_a K-MDP', 'k-means++ K-MDP');
            plot_name_fig = strcat('problems/SONA/results/time-solve.fig');
            plot_name_png = strcat('problems/SONA/results/time-solve.png');
            saveas(gcf, plot_name_fig);
            saveas(gcf, plot_name_png);


    %Represent data data
    
S = readmatrix('problems/SONA/results/states.txt');
S = S(:,1:4);


%Explote data
    
min_dens_aba = S(:,1);
max_dens_aba = S(:,2);
min_abun_so = S(:,3);
max_abun_so = S(:,4);
states = [1:819];


[min_dens_aba_u] = unique(min_dens_aba);
[max_dens_aba_u] = unique(max_dens_aba);
[min_abun_so_u] = unique(min_abun_so);
[max_abun_so_u] = unique(max_abun_so);


 
%Plot state space representation as a matrix.
S_matrix = zeros(length(max_abun_so_u), length(min_dens_aba_u));
it = 1;
for x = 1:length(max_abun_so_u)
    for y = 1:length(min_dens_aba_u)
        S_matrix(x,y) = it;
        it = it + 1;  
    end
end

figure;
surface(min_dens_aba_u, max_abun_so_u, S_matrix);
colorbar('Ticks',[100,200,300,400,500,600,700,800],...
         'TickLabels',{'s_{100}','s_{200}','s_{300}','s_{400}','s_{500}', 's_{600}', 's_{700}', 's_{800}'});
xlabel('Density of northern abalone')
ylabel('Abundance of sea Otter')
title('Sea Otter and Northern Abalone state space');
plot_name_fig = strcat('problems/SONA/results/state-representation.fig');
plot_name_png = strcat('problems/SONA/results/state-representation.png');
saveas(gcf, plot_name_fig);
saveas(gcf, plot_name_png);



%Plot reduced state space representation as a matrix for K=5

S_matrix_reduced = zeros(length(max_abun_so_u), length(min_dens_aba_u));
it = 1;
for x = 1:length(max_abun_so_u)
    for y = 1:length(min_dens_aba_u)
        S_matrix_reduced(x,y) = S2K_astar(it,2);
        it = it + 1;  
    end
end

figure;
surface(min_dens_aba_u, max_abun_so_u, S_matrix, S_matrix_reduced);

un_S_matrix_reduced = [1:1:5];
colormap(cool(length(un_S_matrix_reduced)));
colormap parula;
H = colorbar('Ticks',[1,2,3,4,5],...
         'TickLabels',{'s_{K_1}','s_{K_2}','s_{K_3}','s_{K_4}','s_{K_5}'});
set(H,'ytick',un_S_matrix_reduced);
xlabel('Density of northern abalone')
ylabel('Abundance of sea Otter')
title('Sea Otter and Northern Abalone abstract state space');
plot_name_fig = strcat('problems/SONA/results/abstract-state-representation.fig');
plot_name_png = strcat('problems/SONA/results/abstract-state-representation.png');
saveas(gcf, plot_name_fig);
saveas(gcf, plot_name_png);


%Plot reduced model and policy graphs.


%Name actions
actions = strings(size(P,3), 1);
actions(1) = "I";
actions(2) = "AP";
actions(3) = "C";
actions(4) = "H";
min_prob = 0.001;

%TODO: Check why is crashing here

%[G_reduced_model, p_reduced_model] = model_representation(PK_astar, RK_astar, actions, min_prob);
%[G_reduced_policy, p_reduced_policy] = policy_representation(PK_astar, RK_astar, PolicyK_astar, actions, 0.01);


