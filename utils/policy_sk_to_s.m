function [PolKs] = policy_sk_to_s(K2S, PolicyK, K, NS)
%Converts the PolicyK for K abstracted states to a Policy K for NS original states 


 PolKs=zeros(NS,1);
 for k=1:K
    PolKs(K2S{k})=PolicyK(k);
 end

end