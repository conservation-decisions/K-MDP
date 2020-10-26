%precision
p = 0.000000000000001;
%Discount factor
discount = 0.96;

%load('problems/reserve/P.mat');
%load('problems/reserve/R.mat');
load('P.mat');
load('R.mat');

NS = size(P,1);
NA = size(P,3);

tic;
%Perform value iteration
[Pol]=mdp_value_iteration(P, R, discount);%Get a policy using value iteration
%Evaluate the policy
[V,Q]= mdp_eval_policy_iterative_q(P, R, discount, Pol); % check that this is doing what it is supposed to; Get the value and Q values

time = toc;

K = [20:20:NS];
K = flip(K);
K = [NS, K, 9];

K = [9];%35,30


S = [];

for s1=0:2
    for s2=0:2
        for s3=0:2
            for s4=0:2
                for s5=0:2
                    for s6=0:2
                        for s7=0:2
                            
                           S = [S; s1 s2 s3 s4 s5 s6 s7];
                            
                        end               
                    end
                end
            end
        end
    end
end


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

for i = 1:length(K)
    
    k = K(i);
    
 %fprintf('TRANSITIVE d\n');
 %[PK_Qd, RK_Qd, S2K_Qd{i}, K2S_Qd{i}, PolicyK_Qd, PolKs_Qd, err_transitive(i), t_transitive(i), t_KMDP_transitive(i)] = QdKMDP(k, p, P, R, discount, Q, V);
                
 %fprintf('TRANSITIVE ASTAR\n');
 [PK_astar, RK_astar, S2K_astar, K2S_astar, PolicyK_astar, PolKs_astar, err_astar(i), t_astar(i), t_KMDP_astar(i)] = aStarKMDP(k, p, P, R, discount, V, Pol); 


 %fprintf('MODEL SIMILARITY\n');
 %[PK_m, RK_m, S2K_m, K2S_m, PolicyK_m, PolKs_m, err_m(i), t_m(i), t_KMDP_m(i)] = modelSimilarityKMDP(30, p, P, R, discount, 'yes',V);
    
    
end

K

disp('Q_d')
err_transitive
t_transitive
t_KMDP_transitive

disp('astar')
err_astar
t_astar
t_KMDP_astar

%Store relevant data

save('problems/reserve/results/S2K_Qd.mat', 'S2K_Qd');
save('problems/reserve/results/K2S_Qd.mat', 'K2S_Qd');
save('problems/reserve/results/S2K_astar.mat', 'S2K_astar');
save('problems/reserve/results/K2S_astar.mat', 'K2S_astar');


err_transitive = err_transitive * 100;
err_astar = err_astar * 100;


            figure;
            %plot(K, err_transitive, 'r-o', 'LineWidth', 1)
            hold on;
            plot(K, err_astar, 'b-x', 'LineWidth',1)
            hold on;
            hold off;
            xlabel('K')
            ylabel('gap(%)');
            title('Dynamic reserve site selection');
            legend('Q^*_d K-MDP', 'Q^*_a K-MDP');
            plot_name_fig = strcat('problems/reserve/results/gap_reserve_.fig');
            plot_name_png = strcat('problems/reserve/results/gap_reserve_.png');
            saveas(gcf, plot_name_fig);
            saveas(gcf, plot_name_png);
            
            
            figure;
            %plot(S2K_astar(:,2), S2K_astar(:,1), 'o', 'LineWidth', 1)
            xlabel('S')
            ylabel('S_K');
            title('State abrastraction for K=40');
            %legend('Q^*_d K-MDP', 'Q^*_a K-MDP');
            plot_name_fig = strcat('problems/reserve/results/state_abstraction_.fig');
            plot_name_png = strcat('problems/reserve/results/state_abstraction_.png');
            saveas(gcf, plot_name_fig);
            saveas(gcf, plot_name_png);
            
            
actions = strings(size(P,3), 1);
actions(1) = "R_1";
actions(2) = "R_2";
actions(3) = "R_3";
actions(4) = "R_4";
actions(5) = "R_5";
actions(6) = "R_6";
actions(7) = "R_7";
min_prob = 0.001;
            
            
         %[G_reduced_model, p_reduced_model] = model_representation(PK_astar, RK_astar, actions, min_prob);
         %[G_reduced_policy, p_reduced_policy] = policy_representation(PK_astar, RK_astar, PolicyK_astar, actions, 0.08);
            
