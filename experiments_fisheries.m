%precision
p = 0.00000000000001;
%Discount factor
discount = 0.96;

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

K = [15:10:NS];
K = flip(K);
K = [NS, K, 13];

K = [13];

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


for i = 1:length(K)
    
    k = K(i);
    
    k
    
 fprintf('TRANSITIVE d\n');
 [PK_Qd, RK_Qd, S2K_Qd{i}, K2S_Qd{i}, PolicyK_Qd, PolKs_Qd, err_transitive(i), t_transitive(i), t_KMDP_transitive(i)] = QdKMDP(k, p, P, R, discount, Q, V);
                
 fprintf('TRANSITIVE ASTAR\n');
 [PK_astar, RK_astar, S2K_astar{i}, K2S_astar{i}, PolicyK_astar, PolKs_astar, err_astar(i), t_astar(i), t_KMDP_astar(i)] = aStarKMDP(k, p, P, R, discount, V, Pol); 


 %fprintf('MODEL SIMILARITY\n');
 %[PK_m, RK_m, S2K_m{i}, K2S_m{i}, PolicyK_m, PolKs_m, err_m(i), t_m(i), t_KMDP_m(i)] = modelSimilarityKMDP(k, p, P, R, discount, 'yes',V);
    
    
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

save('problems/fisheries/results/S2K_Qd.mat', 'S2K_Qd');
save('problems/fisheries/results/K2S_Qd.mat', 'K2S_Qd');
save('problems/fisheries/results/S2K_astar.mat', 'S2K_astar');
save('problems/fisheries/results/K2S_astar.mat', 'K2S_astar');


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
            title('Rebuilding global fisheries');
            legend('Q^*_d K-MDP', 'Q^*_a K-MDP');
            plot_name_fig = strcat('problems/fisheries/results/fisheries_gap.fig');
            plot_name_png = strcat('problems/fisheries/results/fisheries_gap.png');
            saveas(gcf, plot_name_fig);
            saveas(gcf, plot_name_png);
            
            
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
            
            
[G_reduced_policy, p_reduced_policy] = policy_representation(PK_astar, RK_astar, PolicyK_astar, actions, 0.08);



