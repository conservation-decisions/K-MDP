K-Markov Decision Process (K-MDP)
Copyright (c) 2020 by Commonwealth Scientific and Industry Research Organization (CSIRO)

Author: Jonathan Ferrer Mestres

The K-MDP toolbox includes a set of algorithms for solving K-MDPs: given an original MDP and a constraint on the 
number of states (K),  generate a reduced state space MDP that  minimizes  the  difference  between  the  original 
optimalMDP value function and the reduced optimal K-MDP value function.

The K-MDP toolbox was developed by the Conservation Decisions team of the Land and Water unit of CSIRO (Australia).

For bug reports and suggestions, please email: jonathan.ferrermestres@csiro.au

[1] Ferrer-Mestres, J.; Dietterich, T. G.; Buffet, O.; and Chades, I.  2020.  Solving K-MDPs. In Proc. International Conference of Automated Planning and Scheduling, 2020.

=================================================
TABLE OF CONTENTS
================================================= 

1. Requirements
2. How to install the K-MDP package
3. How to use the K-MDP package
4. Case studies
5. TODO
6. Acknowledgements


=================================================
REQUIREMENTS
=================================================


MATLAB 
MDPtoolbox

Operating systems:	Ubuntu
			Windows (unofficial)
			Mac OS X (unoficial, not tested)

Tested with MATLAB R2020a and Ubuntu 18.04.



=================================================
HOW TO INSTALL THE K-MDP PACKAGE
=================================================


The K-MDP package requires MATLAB and MDPtoolbox (future updates will make it MDP solver-independent). The current release has been tested with MATLAB R2019a, R2019b and R2020a.
Copy (unziping all the files) or download the content into your system. Make sure you maintain the directory hierarchy. The main directory, K-MDP, should contain 5 subdirectories:

/algorithms	Contains the code of the K-MDP algorithms and the code to build a K-MDP
/kmeans		kmeans++ matlab implementation
/utils		util functions to build K-MDPs
/problems	MDP models defined as reward and transition matrices. Results are saved in this directory.
/MDP_graphs	Code to plot MDP models and policies

Files named experiment_<problem>.m are examples of how to run the K-MDP algorithms and how to plot statistics (gap and times) for different case studies.


=================================================
HOW TO USE THE K-MDP PACKAGE
=================================================

Start MATLAB on the K-MDP directory and set MATLAB path to include K-MDP and all subdirectories. 
MDP models from problems folder are specified as transition matrices (SxSxA) and reward matrices (SxSxA) or (SxA). MDP models are stored as .mat files. If you want to test the provided case studies simply use the experiments_<case_study>.m files or load the .mat files from the case studies folders. If you want to apply the K-MDP algorithms to your MDP models then you have to provide the transition and reward matrices with the previously specified format.

You can choose 3 different algorithms (for more details on how those algorithms work please refer to [1]): Q*_d K-MDP, a*_d K-MDP and k-means K-MDP. All algorithms require solving the MDP beforehand. We use MDPtoolbox to solve the MDPs. Then simply call any of the K-MDP algorithms:

Q*d K-MDP

[AbstractTransitions, AbstractRewards, StoSk, SktoS, AbstractPolicy, AbstractPolicyForOriginalMDP, gap, timeComputingKDMP, timeSolvingKMDP] = QdKMDP(K, precision, Transition, Reward, discount, Q, V);
 
a*d K-MDP

[AbstractTransitions, AbstractRewards, StoSk, SktoS, AbstractPolicy, AbstractPolicyForOriginalMDP, gap, timeComputingKDMP, timeSolvingKMDP] = aStarKMDP(K, precision, Transition, Reward, discount, V, Policy);
 
kmeans K-MDP

[AbstractTransitions, AbstractRewards, StoSk, SktoS, AbstractPolicy, AbstractPolicyForOriginalMDP, gap, timeComputingKDMP, timeSolvingKMDP] = kmeansKMDP(K, Transition, Reward, discount, Q, V);



=================================================
TODO
=================================================

- MDP solver - independent.
- Add the additional K-MDP algorithms into another branch.
- Test it for MAC OS X.
- Test it with Octave.
- Add option to avoid gap computation (if the user only wants the K-MDP, not to compute statistics).

=================================================
ACKNOWLEDGEMENTS
=================================================


This is a joint work with Thomas Dietterich, Olivier Buffet and Iadine Chades.
Special thanks to the Conservation Decisions Team in CSIRO.





