delta = 0.98;
% T: expert, action, current state, future state
% Exp1
T=zeros(4,4,2,2); %added by IC - I was getting an error 26/10/2020 1.19pm BNE, AUS

T(1,1,:,:)= [0.99 0.01;
             0.5 0.5];
T(1,2,:,:)= [0.2 0.8;
             0.05 0.95];      
T(1,3,:,:)= [0.4 0.6;
             0.4 0.6];     
T(1,4,:,:)= [0.6 0.4;
             0.4 0.6];             
% Exp2
T(2,1,:,:)= [0.9 0.1;
             0.6 0.4];
T(2,2,:,:)= [0.55 0.45;
             0.3 0.7];      
T(2,3,:,:)= [0.2 0.8;
             0.1 0.9];     
T(2,4,:,:)= [0.4 0.6;
             0.25 0.75]; 
% Exp3
T(3,1,:,:)= [0.8 0.2;
             0.5 0.5];
T(3,2,:,:)= [0.7 0.3;
             0.4 0.6];      
T(3,3,:,:)= [0.7 0.3;
             0.3 0.7];     
T(3,4,:,:)= [0.5 0.5;
             0.2 0.8]; 
% Exp4
T(4,1,:,:)= [0.95 0.05;
             0.5 0.5];
T(4,2,:,:)= [0.3 0.7;
             0.1 0.9];      
T(4,3,:,:)= [0.5 0.5;
             0.3 0.7];     
T(4,4,:,:)= [0.85 0.15;
             0.4 0.6]; 
P=cell(1,4);
for i=1:4
  P{i} = permute(T(i,:,:,:),[4 3 2 1]);
  P{i}=reshape(P{i},2,8);
end

S=[0;1]; %States
A=[1;2;3;4]; %Actions
X=rectgrid(A,S); %JFM: Actions per state
Ix = getI(X,2);
R=X(:,2);
%inc = 100;
inc = 1;


%Pb    : belief state transition matrix
%   Rb    : belief state reward matrix
%   Sb    : matrix of augmented state variable values



[b,Pb,Rb,Sb,Xb,Ixb]=amdp(inc,P,R,S,X,Ix);
model = struct('P',Pb,'R',Rb,'X',Xb,'Ix',Ixb,'d',delta);
moptions=struct('print',2,'algorithm','i');
results = mdpsolve(model,moptions);
Xopt = results.Xopt;


%%
figure(1); clf
set(gcf,'units','normalized','position',[0.002 0.235 0.983 0.565])

colormap([0.8;0.6;0.4;0.2]+[0 0 0] )
bb=unique(b(:,3));
bb=bb(1:10:101);
k=length(bb);
for i=1:k-1
  subplot(2,k-1,i)
  ind = Xopt(:,2)==0 & Xopt(:,5)==bb(i);
  patchplot(Xopt(ind,3),Xopt(ind,4),Xopt(ind,1),1:4)
  axis square
  if i==1, text(-0.75,0.5,'S=0'); end
  title(['b_3=' num2str(bb(i))])
  xlim([0 1])
  ylim([0 1])
end


for i=1:k-1
  subplot(2,k-1,i+k-1)
  ind = Xopt(:,2)==1 & Xopt(:,5)==bb(i);
  patchplot(Xopt(ind,3),Xopt(ind,4),Xopt(ind,1),1:4)
  axis square
  if i==1, text(-0.75,0.5,'S=1'); end
  title(['b_3=' num2str(bb(i))])
  xlim([0 1])
  ylim([0 1])
end

h=legend('A=1','A=2','A=3','A=4','location','eastoutside');
pos=get(h,'position');
pos(1)=1-1.4*pos(3); pos(2)=0.5-pos(4)/2;
set(h,'position',pos);



           