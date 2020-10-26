% author: iadine.chades@csiro.au
% Maps an optimal action to the vector representing the state of the
% system (sites). Provides the benefits of implementing the action if
% the suggested site was protected successfully as the total number of
% species protected (total_species_p) and the species protected only in the
% designated sites (action).
%
% site = [0 , 0, 1 , 2] means the first 2 sites are available (site(1:2)=0), site 3 is
% reserved (site(3)=1) and site 4 is developped (site(4)=2).

function explore_solution_reserve(init_site,policy,M,P,R)


J=size(M,1);
S=3^J;
T=J+1; % time horizon
action=zeros(T,1);
Tsites=zeros(J,T);
Treward=zeros(T+1,1);
current_site=init_site;

for t=1:T
    s=getState(current_site)+1;
    %Tstates(t,:)=s;
    Tsites(:,t)=current_site';
    action(t)=policy(s);
    if t~= 1
		Treward(t+1)=R(s,action(t))+Treward(t);
    else
		Treward(t+1)=R(s,action(t));
	end
    % Simulating next state s_new and reward associated
    p_s_new = rand(1);
    p = 0;
    s_new = 0;
    while ((p < p_s_new) & (s_new < S))
        s_new = s_new+1;
        p = p + P(s,s_new,action(t));   % be careful. P(s=1,..) corresponds to s=0;
    end;
 
    current_site=getSite(s_new-1,J);
end


figure('color','white');
%subplot(2,1,1);
image(Tsites+1);
cmap=[1, 1, 1
    0.5, 0.5, 0.5;
    0, 0, 0;];
colormap(cmap); 
hcb=colorbar('location','EastOutside','Yticklabel',{'','Available','','Reserved','','Developped',''});
set(hcb,'YTickMode','manual');
xlabel('Time horizon');
ylabel('Sites');
hold on
%subplot(2,1,2);
figure('color','white');
plot(Treward,'color','k');
box off
ylabel('Number of species protected');
xlabel('Time horizon');
end

