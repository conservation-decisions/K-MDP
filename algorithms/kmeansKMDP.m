function [PK, RK, S2K, K2S, PolicyK, PolKs, error, time, time_KMDP] = kmeansKMDP(K, P, R, discount, Q, V)


%kmeansKMDP State abstraction \phi_{a^*_d} K-MDP algorithm
% Arguments -----------------------------------------------
% Let S = number of states, A = number of actions, SK = number of abstract
% states
%   K           = size of the abstracted state space
%   precision   = precision parameter for the binary search
%   P(SxSxA)    = transition matrix 
%              
%   R(SxSxA) or (SxA) = reward matrix
%
%   discount    = discount rate in ]0; 1]
%                 beware to check conditions of convergence for discount = 1.
%   Q(S,A)      = Q values
%   V(S)        = value function
%   metric      = metric for the k-means++ clustering technique. 
%
% Output ----------------------------------------------------
%   PK(SKxSKxA) = transition matrix for the K-MDP model.
%   RK(SkxSKxA) or (SKxA)= reward matrix for the K-MDP model. 
%   S2K(S)      = mapping from state S to its abstract state SK
%   K2S(SK)     = mapping fro abstract state SK to state S
%   PolicyK(SK) = abstract optimal policy
%   PolKs(S)    = abstract optimal policy to state space S
%   error       = Gap between the original optimal value function and the
%                 abstract optimal value function
%   time        = Time required to compute the K state abstraction
%   time_KMDP   = Time required to solve the K-MDP model using value
%                 iteration


    %Number of states NS
    NS = size(P,1);
                    
    tic;
    

    
    [L,C]=kmeanspp(Q,K); % change that to account for Q* metric; Abstracted state space
    


    hold on;

    [PK,RK,S2K,K2S]=buildKMDP(P,R,L,K); % done might need to normalise; Build the KMDP with k states and state space L              
                
    time = toc;
                     
    tic;

    % solve the KMDP
    [PolicyK]=mdp_value_iteration(PK,RK,discount);% Get the policy for the abstracted KMDP               

    %Converts the PolicyK for K abstracted states to a Policy K for NS original states
    [PolKs] = policy_sk_to_s(K2S, PolicyK, K, NS);
                 
    %Value of applying PolicyK to the original MDP
    VK = mdp_eval_policy_iterative(P, R, discount, PolKs);                
                
    time_KMDP = toc;
                                
    %Compute the error between both policies        
                
    diff_policies = V - VK;  
    [error, index] = max(diff_policies);              
    error = error/V(index);
           


end
    