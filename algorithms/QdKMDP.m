function [PK, RK, S2K, K2S, PolicyK, PolKs, error, time, time_KMDP] = QdKMDP(K, precision, P, R, discount, Q, V)

% QdKMDP State abstraction \phi_{a^*_d} K-MDP algorithm
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
%   V(S)        = value function
%   Policy(S)   = optimal policy
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
    %Number of actions
    NA = size(P,3);

    VMAX = max(V); %Take the maximum value
    
    tic;
    
    d_lower = 0;
    d_upper = VMAX;

    p = d_upper - d_lower; %To refine precision 

    %Create a structure to store the bindings
    bindings = zeros(NS, 2);

    while p > precision


          p = d_upper - d_lower; %To refine precision  
          d = d_lower + (d_upper - d_lower)/2; %Binary search on d


          %Compute the ceil of all Q values
          for s=1:NS
              for a=1:NA
                 bindings(s,a) = ceil(Q(s,a)/d);   
              end 
          end        


          [B,ia,ic] = unique(bindings, 'rows');


          tmp_number_abstractions = length(ia);


          if tmp_number_abstractions <= K
              %Here K must be the number of clusters in L. In this case is the size of ic;
              L = ic';
              number_abstractions = tmp_number_abstractions;
              d_upper = d;
          else
              d_lower = d;
          end

    end
    
    if number_abstractions > K
        
        PK = NaN;
        RK = NaN;
        S2K = NaN;
        K2S = NaN;
        PolicyK = NaN;
        PolKs = NaN;
        error = NaN;
        time = NaN;
        time_KMDP = NaN;
        
        return;
        
    end
    
    

    %Build the KMDP with the abstracted state space
   [PK,RK,S2K,K2S]=buildKMDP(P, R, L, number_abstractions); % done might need to normalise; Build the KMDP with k states and state space L    

    time = toc;
    
    tic;
          
    % Solve the KMDP
    [PolicyK]=mdp_value_iteration(PK,RK,discount);% Get the policy for the abstracted KMDP

    %Converts the PolicyK for K abstracted states to a Policy K for NS original states
    [PolKs] = policy_sk_to_s(K2S, PolicyK, number_abstractions, NS);

    %Value of applying PolicyK to the original MDP
    VK = mdp_eval_policy_iterative(P, R, discount, PolKs);
     
    time_KMDP = toc;
         
    %Compute the error between both policies
    diff_policies = V - VK;
      
    [error, index] = max(diff_policies);
                      
    error = error/V(index);
    
            
    if V(index) == 0
        error = 0;
    end
                 
    
 

end