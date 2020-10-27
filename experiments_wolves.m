%Add dependencies
addpath('kmeans')
addpath('utils')
addpath('algorithms')
addpath('MDP_graph_representation')

%precision
p = 0.00000000001;
%Discount factor
discount = 0.96;

min_pop = 50;
max_pop = 250;

NS = 250;
NA = 251;


%Not available yet given the size of the mat files.

P_mat = strcat("problems/wolves/wolves_P_",num2str(NS),"_",num2str(NA),"_min_",num2str(min_pop),"_max_",num2str(max_pop),".mat");
R_mat = strcat("problems/wolves/wolves_R_",num2str(NS),"_",num2str(NA),"_min_",num2str(min_pop),"_max_",num2str(max_pop),".mat");


load(P_mat);
load(R_mat);

NS = size(P,1);

%Solve the MDP

tic;

%Perform value iteration
[Pol]=mdp_value_iteration(P, R, discount);%Get a policy using value iteration

%Evaluate the policy
[V,Q]= mdp_eval_policy_iterative_q(P, R, discount, Pol); % check that this is doing what it is supposed to; Get the value and Q values

time_MDP = toc;


K = [10:10:NS];


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

save('problems/wolves/results/S2K_Qd.mat', 'S2K_Qd');
save('problems/wolves/results/K2S_Qd.mat', 'K2S_Qd');
save('problems/wolves/results/S2K_astar.mat', 'S2K_astar');
save('problems/wolves/results/K2S_astar.mat', 'K2S_astar');
save('problems/wolves/results/S2K_kmeans.mat', 'S2K_kmeans');
save('problems/wolves/results/K2S_kmeans.mat', 'K2S_kmeans');


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
            title('Wolf culling error');
            legend('Q^*_d K-MDP', 'Q^*_a K-MDP', 'k-means++ K-MDP');
            plot_name_fig = strcat('problems/wolves/results/wolves_gap_250_251.fig');
            plot_name_png = strcat('problems/wolves/results/wolves_gap_250_251.png');
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
            title('Wolf culling time to compute the K-MDP');
              legend('Q^*_d K-MDP', 'Q^*_a K-MDP', 'k-means++ K-MDP');
            plot_name_fig = strcat('problems/wolves/results/wolves_time-compute_250_251.fig');
            plot_name_png = strcat('problems/wolves/results/wolves_time-compute_250_251.png');
            saveas(gcf, plot_name_fig);
            saveas(gcf, plot_name_png);
            
            
            %Time to solveKDMP
            
            figure;
            plot(K, t_KMDP_transitive, 'r-o', 'LineWidth', 1)
            hold on;
            plot(K, t_KMDP_astar, 'b-x', 'LineWidth',1)
            hold on;
            plot(K, t_KMDP_kmeans, 'b-x', 'LineWidth',1)
            hold off;
            xlabel('K')
            ylabel('time(sec.)');
            title('Wolf culling time to solve the K-MDP');
              legend('Q^*_d K-MDP', 'Q^*_a K-MDP', 'k-means++ K-MDP');
            plot_name_fig = strcat('problems/wolves/results/wolves_time-solve_250_251.fig');
            plot_name_png = strcat('problems/wolves/results/wolves_time-solve_250_251.png');
            saveas(gcf, plot_name_fig);
            saveas(gcf, plot_name_png);
            
            
            %Plot reduced model and policy graphs.

            
            
actions = strings(size(P,3), 1);
actions(1) = "0%_H";
actions(2) = "10%_H";
actions(3) = "20%_H";
actions(4) = "30%_H";
actions(5) = "40%_H";
actions(6) = "50%_H";
actions(7) = "60%_H";
actions(8) = "70%_H";
actions(9) = "80%_H";
actions(10) = "90%_H";
actions(11) = "100%_H";
min_prob = 0.001;

%TODO: Check why is crashing here

            

%[G_reduced_model, p_reduced_model] = model_representation(PK_astar, RK_astar, actions, min_prob);
%[G_reduced_policy, p_reduced_policy] = policy_representation(PK_astar, RK_astar, PolicyK_astar, actions, 0.01);

