function [P, R] = mdp_example_reserve (M, pj)


% mdp_example_reserve   Generate a Markov Decision Process example based on
%                       a simple reserve design problem
%                       (see the related documentation for more detail)
% Arguments -------------------------------------------------------------
%   M(JxI) = species distribution across sites
%        J = number of sites (> 0), optional (default 5)
%        I = number of species (>0), optional (default 7)
%   pj = probability of development occurence, in ]0, 1[, optional (default 0.1)
% Evaluation -------------------------------------------------------------
%   P(SxSxA) = transition probability matrix
%   R(SxA) = reward matrix
%   M(JxI) = random species distribution across sites

% MDPtoolbox: Markov Decision Processes Toolbox
% Copyright (C) 2009  INRA
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%    * Neither the name of the <ORGANIZATION> nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
% OF THE POSSIBILITY OF SUCH DAMAGE.


% arguments checking
if nargin >= 1 & size(M)<= 1
    disp('----------------------------------------------------------')
    disp('MDP Toolbox ERROR: M is a JxI matrix with Number of sites J and species I must be greater than 0')
    disp('----------------------------------------------------------')
elseif nargin >= 2 & ( pj < 0 | pj > 1 )
    disp('-----------------------------------------------------------')
    disp('MDP Toolbox ERROR: Probability pj must be in [0; 1]')
    disp('-----------------------------------------------------------')
else

    % initialization of optional arguments
    if nargin < 2; pj=0.1; end;
    if nargin < 1
        J=5;I=7;
        M=round(rand(J,I)); % currently random but feel free to provide a matrix as an input variableI=5; end;
    end;
    [J,I]=size(M);
    % Definition of states
    S=3^J;    % for each site, 0 is available; 1 is reserved; 2 is developed ;
    A=J;      % action space

    % There are J actions corresponding to the selection of a site for
    % reservation. A site can only be reserved if it is available.
    % By convention we will use a ternary base where state #0 is the state
    % that corresponds to [0,0, ..,0] all sites are available. State #1 is
    % [1,0,0,...,0]; state 2 is [2,0,0 ..,0] and state 3 is [0,1,0, .. 0] and so forth.
    % for example 
    % site = [0,0,1,2] means the first 2 sites are available (site(1:2)=0), site 3 is
    % reserved (site(3)=1) and site 4 is developped (site(4)=2).
    
    % Build P(AxSxS)
    % complexity is in SxAx2^navail; with 2^navail<=S;
    
    P=zeros(S,S,A);
    for s1=1:S  % for all states
        site1= getSite(s1-1,J); % state n becomes n-1 // MATLAB index starts at 1 not 0
        for a=1:A   % for all actions
            site2=site1; % site2 represents the site after action a is performed
            if site1(a)==0  % if a is an available site
                site2(a)=1; % site a is reserved
            end
            availSite=find(site2 == 0); % where are the potential available sites?
            if (~isempty(availSite)) % some sites might become candidate for development
                navail=length(availSite);   % how many?
                siten=ones(2^navail,1)*site2; % siten is the set of successors state, number of successors = 2^ navail sites
                aux=zeros(1,navail);          % trick to build the set of successor states
                for k=1:2^(navail)      % there are exactly 2^navail successors
                    siten(k,availSite)=aux.*2;   % init to aux *.2 because developped site are #2;
                    ndev= sum(abs(site2(availSite)-aux)); % how many are developped
                    s2=getState(siten(k,:))+1;  % corresponding state
                    P(s1,s2,a)=pj^ndev*(1-pj)^(navail-ndev);    % corresponding prob of transition
                    aux=dec2binvec(binvec2dec(aux)+1, navail);  % aux <- aux + 1 (binary addition)!
                                        
                end
            else % only one future state possible 
               s2=getState(site2)+1;
               P(s1,s2,a)=1; 
            end
        end
    end  
    
    % build R
    
    R=zeros(S,A);
    for s1=1:S
        site1= getSite(s1-1,J); % state n becomes n-1 // MATLAB index starts at 1 not 0
        for a=1:J
            if (site1(a)~=2)
                targetSp=M(a,:);
                reservedSites=find(site1==1);
                if ~isempty(reservedSites)
                    reservedSp=max(M(reservedSites,:),[],1);
                    R(s1,a)=sum(max(0.0001,targetSp-reservedSp));
                else
                    R(s1,a)=sum(max(0.0001,targetSp));
                end
            end
        end
    end
end