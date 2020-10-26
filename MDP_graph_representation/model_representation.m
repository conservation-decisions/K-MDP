function [G,p] = policy_representation(P, R, actions, min_prob, rewards)

    

    %Change Nan for 0s
    P(isnan(P)) = 0;
    
    NS = size(P,1);
    NA = size(P, 3);
    
    
    if nargin < 3
      actions = strings(NA, 1);
          for a = 1:NA
              actions(a) = num2str(a);
          end
      min_prob = 0.0;
      rewards = 'n';
    end
    
    if nargin < 4
        min_prob = 0.0;
        rewards = 'n';
    end
    
     if nargin < 5
        rewards = 'n';
     end
   
   
     
    init = [];
    goal = [];

    elabels = {};
        
    for s = 1:NS
                    
        for sp = 1:NS
           
            for a = 1:NA
            
                if P(s,sp,a) >= min_prob
                    init = [init s];
                    goal = [goal sp];

                    elabels = [elabels strcat(actions(a), '(', num2str(round(P(s,sp,a),3)), ')')];
                end
            
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
     
    plot_name_fig = strcat('model_', num2str(size(P,1)) ,'.fig');
    plot_name_png = strcat('model_', num2str(size(P,1)) ,'.png');
    saveas(gcf, plot_name_fig);
    saveas(gcf, plot_name_png);
     
     
       
end