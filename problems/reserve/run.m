% Add the necessary paths
%addpath('')

J = 7;
I = 20;


% Generate the species richness matrix
rand('seed',0)
M=round(rand(J,I));

discount = 0.96;

% Generate the transition and reward matrix
[P, R] = mdp_example_reserve(M, 0.2);

R = R - min(R(:));
R = R ./ max(R(:));


save('P.mat', 'P');
save('R.mat', 'R');


% Solve the reserve design problem
[policy, iter, cpu_time] = mdp_value_iteration(P, R, discount, 0.001);

[V,Q]= mdp_eval_policy_iterative_q(P, R, discount, policy); % check that this is doing what it is supposed to; Get the value and Q values


% Explore solution with initial state all sites available
%explore_solution_reserve([0 0 0 0 0 0 0],policy,M,P,R)


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
