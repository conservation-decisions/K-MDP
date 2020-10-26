function [G,p] = policy_representation(P, R, Policy, actions, min_prob, rewards, prune)

    
    %Change Nan for 0s
    P(isnan(P)) = 0;
    
    NS = size(P,1);
    NA = size(P, 3);
    
    Policy_size = length(Policy);
    
      if nargin < 4
        actions = strings(NA, 1);
          for a = 1:NA
              actions(a) = num2str(a);
          end
        min_prob = 0.0;
        rewards = 'n';
        prune = 'n';
      end
    
    if nargin < 5
        min_prob = 0.0;
        rewards = 'n';
        prune = 'n';
    end
    
     if nargin < 6
        rewards = 'n';
        prune = 'n';
     end
     
     if nargin < 7
        prune = 'n';
     end
   
   
     
    init = [];
    goal = [];

    elabels = {};
        
    for s = 1:length(Policy)
        
        action = Policy(s);
       
        for sp = 1:NS
           
            
            if P(s,sp,action) >= min_prob 
                init = [init s];
                goal = [goal sp];
                
                elabels = [elabels strcat(actions(action), '(', num2str(round(P(s,sp,action),3)), ')')];
            end
            
        end
        
    end
     
    
     G = digraph(init, goal);
      
    nodes = strings(NS,1);
    for n = 1:NS
        
        nodes(n) = strcat('S_K_{',num2str(n),'}');
    end
     
    
    G.Nodes.Name = nodes;
    
    
    p = plot(G,'NodeLabel',nodes,'EdgeLabel',elabels','LineWidth',2, 'ArrowSize', 15)
    p.NodeColor = 'r';
    p.MarkerSize = 25;
    p.NodeFontWeight = 'bold';
    p.NodeFontSize = 12;
    p.EdgeFontWeight = 'bold';
    p.EdgeFontSize = 9;
    p.LineWidth = 2;
    
     
    plot_name_fig = strcat('policy_', num2str(size(P,1)) ,'.fig');
    plot_name_png = strcat('policy_', num2str(size(P,1)) ,'.png');
    saveas(gcf, plot_name_fig);
    saveas(gcf, plot_name_png);
     
     
       

end