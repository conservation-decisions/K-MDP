%precision
p = 0.0000000000001;
%Discount factor
discount = 0.96;

min_pop = 50;
max_pop = 250;

NS = 250;
NA = 251;

%P_mat = strcat("problems/wolves/min_",num2str(min_pop),"_max_",num2str(max_pop),"/wolves_P_",num2str(NS),"_",num2str(NA),"_min_",num2str(min_pop),"_max_",num2str(max_pop),".mat");
%R_mat = strcat("problems/wolves/min_",num2str(min_pop),"_max_",num2str(max_pop),"/wolves_R_",num2str(NS),"_",num2str(NA),"_min_",num2str(min_pop),"_max_",num2str(max_pop),".mat");

P_mat = strcat("problems/wolves/wolves_P_",num2str(NS),"_",num2str(NA),"_min_",num2str(min_pop),"_max_",num2str(max_pop),".mat");
R_mat = strcat("problems/wolves/wolves_R_",num2str(NS),"_",num2str(NA),"_min_",num2str(min_pop),"_max_",num2str(max_pop),".mat");


%load('problems/wolves/wolves_P_251_11_min_5_max_15.mat');
%load('problems/wolves/wolves_R_251_11_min_5_max_15.mat');
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

%K
%K = [200:200:NS];
%K = flip(K);
%K = [K, 6];
%K = [50:5:NS];
%K = flip(K);
%K = [K, 6];
K = [50];

%K = [50];

%Initialize structures for statistics extraction
err_transitive = zeros(length(K),1);
t_transitive = zeros(length(K),1);
t_KMDP_transitive = zeros(length(K),1);

err_astar = zeros(length(K),1);
t_astar = zeros(length(K),1);
t_KMDP_astar = zeros(length(K),1);

err_m = zeros(length(K),1);
t_m = zeros(length(K),1);
t_KMDP_m = zeros(length(K),1);

S2K_Qd = cell(length(K),1);
K2S_Qd = cell(length(K),1);

S2K_astar = cell(length(K),1);
K2S_astar = cell(length(K),1);

S2K_m = cell(length(K),1);
K2S_m = cell(length(K),1);


%[PK_Qd, RK_Qd, S2K, K2S, PolicyK_Qd, PolKs_Qd, err_t, t_t, t_KMDP_t] = QdKMDP(K, p, P, R, discount, Q, V);
%[PK_astar, RK_astar, S2K_a, K2S_a, PolicyK_a, PolKs_a, err_a, t_a, t_KMDP_a] = aStarKMDP(K, p, P, R, discount, V, Pol); 

%[PK_m, RK_m, S2K_m, K2S_m, PolicyK_m, PolKs_m, err_m, t_m, t_KMDP_m] = modelSimilarityKMDP(847, p, P, R, discount, 'yes',V);
%clear V;


%Perform value iteration
%[PolicyKM]=mdp_value_iteration(PK_m, RK_m, discount);%Get a policy using value iteration

%Evaluate the policy
%[VKM,QKM]= mdp_eval_policy_iterative_q(PK_m, RK_m, discount, PolicyKM); % check that this is doing what it is supposed to; Get the value and Q values


%[PK_Qd, RK_Qd, S2K_Qd, K2S_Qd, PolicyK_Qd, PolKs_Qd, err_transitive, t_transitive, t_KMDP_transitive] = QdKMDP(5, p, PK_m, RK_m, discount, QKM, VKM);
%[PK_astar, RK_astar, S2K_astar, K2S_astar, PolicyK_astar, PolKs_astar, err_astar, t_astar, t_KMDP_astar] = aStarKMDP(50, p, PK_m, RK_m, discount, VKM, PolicyKM); 




for i = 1:length(K)
    
    k = K(i);
    
 fprintf('TRANSITIVE d\n');
% [PK_Qd, RK_Qd, S2K_Qd{i}, K2S_Qd{i}, PolicyK_Qd, PolKs_Qd, err_transitive(i), t_transitive(i), t_KMDP_transitive(i)] = QdKMDP(k, p, P, R, discount, Q, V);
                
 fprintf('TRANSITIVE ASTAR\n');
 [PK_astar, RK_astar, S2K_astar{i}, K2S_astar{i}, PolicyK_astar, PolKs_astar, err_astar(i), t_astar(i), t_KMDP_astar(i)] = aStarKMDP(k, p, P, R, discount, V, Pol); 


 %fprintf('MODEL SIMILARITY\n');
 %[PK_m, RK_m, S2K_m{i}, K2S_m{i}, PolicyK_m, PolKs_m, err_m(i), t_m(i), t_KMDP_m(i)] = modelSimilarityKMDP(817, p, P, R, discount, 'yes',V);
    
    
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

%Store relevant data

save('problems/wolves/results/S2K_Qd_250_251.mat', 'S2K_Qd');
save('problems/wolves/results/K2S_Qd_250_251.mat', 'K2S_Qd');
save('problems/wolves/results/S2K_astar_250_251.mat', 'S2K_astar');
save('problems/wolves/results/K2S_astar_250_251.mat', 'K2S_astar');


err_transitive = err_transitive * 100;
err_astar = err_astar * 100;


%Plot releveant data


%Error

            figure;
            plot(K, err_transitive, 'r-o', 'LineWidth', 1)
            hold on;
            plot(K, err_astar, 'b-x', 'LineWidth',1)
            hold on;
            hold off;
            xlabel('K')
            ylabel('gap(%)');
            title('Wolf culling error');
            legend('Q^*_d K-MDP', 'Q^*_a K-MDP');
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
            hold off;
            xlabel('K')
            ylabel('time(sec.)');
            title('Wolf culling time to compute the K-MDP');
            legend('Q^*_d K-MDP', 'Q^*_a K-MDP');
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
            hold off;
            xlabel('K')
            ylabel('time(sec.)');
            title('Wolf culling time to solve the K-MDP');
            legend('Q^*_d K-MDP', 'Q^*_a K-MDP');
            plot_name_fig = strcat('problems/wolves/results/wolves_time-solve_250_251.fig');
            plot_name_png = strcat('problems/wolves/results/wolves_time-solve_250_251.png');
            saveas(gcf, plot_name_fig);
            saveas(gcf, plot_name_png);
            
            
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
            

%[G_reduced_model, p_reduced_model] = model_representation(PK_astar, RK_astar, actions, min_prob);
%[G_reduced_policy, p_reduced_policy] = policy_representation(PK_astar, RK_astar, PolicyK_astar, actions, 0.01);

