%Add dependencies
addpath('kmeans')
addpath('utils')
addpath('algorithms')
addpath('MDP_graph_representation')


%precision
p = 0.00000000000001;
%Discount factor
discount = 0.96;


%Choose the MDP model from the ones available on problems/fisheries.
NS = 1001;
NA = 11;

P_mat = strcat("problems/fisheries/fisheries_P_",num2str(NS),"_",num2str(NA),".mat");
R_mat = strcat("problems/fisheries/fisheries_R_",num2str(NS),"_",num2str(NA),".mat");

load(P_mat);
load(R_mat);


%Solve the MDP

tic;

%Perform value iteration
[Pol]=mdp_value_iteration(P, R, discount);%Get a policy using value iteration

%Evaluate the policy
[V,Q]= mdp_eval_policy_iterative_q(P, R, discount, Pol); % check that this is doing what it is supposed to; Get the value and Q values

Q(isnan(Q))=0;

time_MDP = toc;

%K
K = [15:10:NS];
K = flip(K);
K = [NS, K, 13];

K = 13;


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

S2K_kmeans = cell(length(K),1);
K2S_kmeans = cell(length(K),1);


for i = 1:length(K)
    
    k = K(i);
    
 fprintf('Qd\n');
 [PK_Qd, RK_Qd, S2K_Qd{i}, K2S_Qd{i}, PolicyK_Qd, PolKs_Qd, err_transitive(i), t_transitive(i), t_KMDP_transitive(i)] = QdKMDP(k, p, P, R, discount, Q, V);
                
 fprintf('a*_d\n');
 [PK_astar, RK_astar, S2K_astar{i}, K2S_astar{i}, PolicyK_astar, PolKs_astar, err_astar(i), t_astar(i), t_KMDP_astar(i)] = aStarKMDP(k, p, P, R, discount, V, Pol); 

 fprintf('K-means \n');
 [PK_kmeans, RK_Kmeans, S2K_kmeans{i}, K2S_kmeans{i}, PolicyK_kmeans, PolKs_kmeans, err_kmeans(i), t_kmeans(i), t_KMDP_kmeans(i)] = kmeansKMDP(k, P, R, discount, Q, V);
        
    
end


%DATA to plot and extract


K

disp('Q_d')
err_transitive'
t_transitive'
t_KMDP_transitive'

disp('astar')
err_astar'
t_astar'
t_KMDP_astar'

disp('kmeans++')
err_kmeans
t_kmeans
t_KMDP_kmeans


%Store relevant data

save('problems/fisheries/results/S2K_Qd.mat', 'S2K_Qd');
save('problems/fisheries/results/K2S_Qd.mat', 'K2S_Qd');
save('problems/fisheries/results/S2K_astar.mat', 'S2K_astar');
save('problems/fisheries/results/K2S_astar.mat', 'K2S_astar');
save('problems/fisheries/results/S2K_kmeans.mat', 'S2K_kmeans');
save('problems/fisheries/results/K2S_kmeans_250_251.mat', 'K2S_kmeans');


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
            title('Rebuilding global fisheries');
            legend('Q^*_d K-MDP', 'Q^*_a K-MDP', 'k-means++ K-MDP');
            plot_name_fig = strcat('problems/fisheries/results/fisheries_gap.fig');
            plot_name_png = strcat('problems/fisheries/results/fisheries_gap.png');
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
            plot_name_fig = strcat('problems/fisheries/results/time-compute.fig');
            plot_name_png = strcat('problems/fisheries/results/time-compute.png');
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
            plot_name_fig = strcat('problems/fisheries/results/time-solve.fig');
            plot_name_png = strcat('problems/fisheries/results/time-solve.png');
            saveas(gcf, plot_name_fig);
            saveas(gcf, plot_name_png);
            
            
            %Plot reduced model and policy graphs.
            
            
            
actions = strings(size(P,3), 1);
actions(1) = "0_H";
actions(2) = "100_H";
actions(3) = "200_H";
actions(4) = "300_H";
actions(5) = "400_H";
actions(6) = "500_H";
actions(7) = "600_H";
actions(8) = "700_H";
actions(9) = "800_H";
actions(10) = "900_H";
actions(11) = "1000_H";
min_prob = 0.001;
            
            
%TODO: Check why is crashing here

%[G_reduced_model, p_reduced_model] = model_representation(PK_astar, RK_astar, actions, min_prob);
%[G_reduced_policy, p_reduced_policy] = policy_representation(PK_astar, RK_astar, PolicyK_astar, actions, 0.08);



