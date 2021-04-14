% Noncontextuality inequality involving nontrivial preparation and measurement operational equivalences, from
% D. Schmid, R. W. Spekkens, E. Wolfe, Phys. Rev. A 97, 062103 (2018).
% The scenario is outlined in Section VI.C, and here we look at all the inequalities in Eq. (26)

%% Specify the number of preparations, measurements and outcomes per measurements
% The preparations are labelled as in the paper cited above. The POVM elements are labelled naturally as
% 1 -> [0|M_1], 2 -> [1|M_1], 3 -> [0|M_2], etc.
nX = 6;
nY = 3; 
nB = 2; 

%% Specifies the score function for the task, so that score = sum_{b,x,y} w(b,x,y)*p(b|x,y)
% Inequality 1: p(0|1,1) ≤ 1
w1 = zeros(nB,nX,nY);
w1(1,1,1) = 1;
% Inequality 2: p(0|1,1) + p(0|3,2) + p(0|5,3) ≤ 2.5
w2 = zeros(nB,nX,nY);
w2(1,1,1) = 1;
w2(1,3,2) = 1;
w2(1,5,3) = 1;
% Inequality 3: p(0|1,1) + p(0|2,2) + p(0|5,3) ≤ 2.5
w3 = zeros(nB,nX,nY);
w3(1,1,1) = 1;
w3(1,2,2) = 1;
w3(1,5,3) = 1;
% Inequality 4: p(0|1,1) - p(0|4,1) - 2p(0|5,1) - 2p(0|2,2) + 2p(0|3,2) + 2p(0|5,3)  ≤ 3
w4 = zeros(nB,nX,nY);
w4(1,1,1) = 1;
w4(1,4,1) = -1;
w4(1,5,1) = -2;
w4(1,2,2) = -2;
w4(1,3,2) = 2;
w4(1,5,3) = 2;
% Inequality 5: 2p(0|1,1) - p(0|2,2) + 2p(0|3,2) ≤ 3
w5 = zeros(nB,nX,nY);
w5(1,1,1) = 2;
w5(1,2,2) = -1;
w5(1,3,2) = 2;
% Inequality 6: p(0|1,1) - p(0|5,1) + p(0|2,2) + p(0|3,2) + 2p(0|5,3) ≤ 4
w6 = zeros(nB,nX,nY);
w6(1,1,1) = 1;
w6(1,5,1) = -1;
w6(1,2,2) = 1;
w6(1,3,2) = 1;
w6(1,5,3) = 2;
% Inequality 7: p(0|1,1) - p(0|5,1) + 2p(0|2,2) + 2p(0|5,3) ≤ 4
w7 = zeros(nB,nX,nY);
w7(1,1,1) = 1;
w7(1,5,1) = -1;
w7(1,2,2) = 2;
w7(1,5,3) = 2;

NC_bounds = [1, 2.5, 2.5, 3, 3, 4, 4];
Q_bounds = [1, 3, 2.8660, 3.5, 3.3660, 4.6889, 4.6458];

% Constraints on probability distribution: For this inequality we don't have any.
p_constraints = [];

%% Specify the operational equivalences:
% Each operational equivalence is a set (cell) of 2x|S_k| matrices
% The first row specifies the sets (S_k), the second the weights for the corresponding preparations (the \xi_k)
Popequiv = cell(3,1);
Popequiv{1} = [1 2; 1/2 1/2];
Popequiv{2} = [3 4; 1/2 1/2];
Popequiv{3} = [5 6; 1/2 1/2];
E_P = {Popequiv};

% The operational equivalences for the measurement effects are specified analagously
Mopequiv = cell(2,1);
Mopequiv{1} = [1 3 5; 1/3 1/3 1/3];
Mopequiv{2} = [2 4 6; 1/3 1/3 1/3];
E_M = {Mopequiv};

%% Now we specify the levels of the hierarchy we want to use.
% 1: states; 2: measurements, 3: sigma (here we have only one (preparation) operational equivalence)
% E.g., [1 2 2] corresponds to including all rho*M*M, terms in the monomial list \mathcal{S}
% or, in the case of the sqrt(E_{b|m}) variant, it means rho*\sqrt{E_{b|y}}*sqrt{E_{b|y}}, etc.

% Main moment matrix:
% ====== Level 1 (all) ======
levelsS = {1, 2, 3};
% ====== Level 2 terms ======
levelsS{end+1} = [1, 1];
levelsS{end+1} = [1, 2];
levelsS{end+1} = [1, 3];
levelsS{end+1} = [2, 1];
levelsS{end+1} = [2, 2];
levelsS{end+1} = [2, 3];
% ====== Level 3 terms ======
levelsS{end+1} = [1, 1, 1];
levelsS{end+1} = [1, 1, 2];
levelsS{end+1} = [1, 1, 3];
levelsS{end+1} = [2, 2, 1];
levelsS{end+1} = [2, 2, 2];
levelsS{end+1} = [2, 2, 3];

% Localising matrices (L for the operational equivalences, O for the operator positivity):
% ====== Level 1 terms ======
levelsL = {1, 2};
% ====== Level 2 terms ======
levelsL{end+1} = [1, 1];
levelsL{end+1} = [2, 2];

% Use same monomial list for the localising matrices for operator positivity
levelsO = levelsL;

levels = {levelsS, levelsL, levelsO};

%% Specify the symmetries in both the states and the measurements (to be understood as simultaneous with those for the states)
% nX preparations rho_1,...,rho_nX ~= [1,...,nX]
% nY*nB measurements M_{1|1},...,M_{nB|nY} ~= [1,...,nY*nB]

% Inequality 1: p(0|1,1) ≤ 1
stateSymGens1 = cell(3,1);
measurementSymGens1 = cell(3,1);
% P3 <-> P4
stateSymGens1{1} = [1 2 4 3 5 6];
measurementSymGens1{1} = [1 2 3 4 5 6];
% (P3,P4) <-> (P5,P6)
stateSymGens1{2} = [1 2 5 6 3 4];
measurementSymGens1{2} = [1 2 3 4 5 6];
% M2 <-> M3
stateSymGens1{3} = [1 2 3 4 5 6];
measurementSymGens1{3} = [1 2 5 6 3 4];
symmetryGenerators1 = {stateSymGens1, measurementSymGens1};

% Inequality 2: p(0|1,1) + p(0|3,2) + p(0|5,3) ≤ 2.5
stateSymGens2 = cell(2,1);
measurementSymGens2 = cell(2,1);
% p(0|1,1) <-> p(0|3,2), while keeping OEs satisfied
stateSymGens2{1} = [3 4 1 2 5 6];
measurementSymGens2{1} = [3 4 1 2 5 6];
% p(0|1,1) <-> p(0|5,3), while keeping OEs satisfied
stateSymGens2{2} = [5 6 3 4 1 2];
measurementSymGens2{2} = [5 6 3 4 1 2];
symmetryGenerators2 = {stateSymGens2, measurementSymGens2};

% Inequality 3: p(0|1,1) + p(0|2,2) + p(0|5,3) ≤ 2.5
stateSymGens3 = cell(2,1);
measurementSymGens3 = cell(2,1);
% p(0|1,1) <-> p(0|2,2), while keeping OEs satisfied
stateSymGens3{1} = [2 1 3 4 5 6];
measurementSymGens3{1} = [3 4 1 2 5 6];
% P3 <-> P4
stateSymGens3{2} = [1 2 4 3 5 6];
measurementSymGens3{2} = [1 2 3 4 5 6];
symmetryGenerators3 = {stateSymGens3, measurementSymGens3};

% Inequality 4: p(0|1,1) - p(0|4,1) - 2p(0|5,1) - 2p(0|2,2) + 2p(0|3,2) + 2p(0|5,3)  ≤ 3
% No symmetries here that preserve the operational equivalences and objective function :-(
symmetryGenerators4 = {{[1 2 3 4 5 6]}, {[1 2 3 4 5 6]}};

% Inequality 5: 2p(0|1,1) - p(0|2,2) + 2p(0|3,2) ≤ 3
stateSymGens5 = cell(1,1);
measurementSymGens5 = cell(1,1);
% P5 <-> P6
stateSymGens5{1} = [1 2 3 4 6 5];
measurementSymGens5{1} = [1 2 3 4 5 6];
symmetryGenerators5 = {stateSymGens5, measurementSymGens5};

% Inequality 6: p(0|1,1) - p(0|5,1) + p(0|2,2) + p(0|3,2) + 2p(0|5,3) ≤ 4
% No symmetries here that preserve the operational equivalences and objective function :-(
symmetryGenerators6 = {{[1 2 3 4 5 6]}, {[1 2 3 4 5 6]}};

% Inequality 7: p(0|1,1) - p(0|5,1) + 2p(0|2,2) + 2p(0|5,3) ≤ 4
stateSymGens7 = cell(1,1);
measurementSymGens7 = cell(1,1);
% P3 <-> P4
stateSymGens7{1} = [1 2 4 3 5 6];
measurementSymGens7{1} = [1 2 3 4 5 6];
symmetryGenerators7 = {stateSymGens7, measurementSymGens7};

%% Setup options for hierarchy and solve
options = struct();
options.verbose = true;

options.pure = false; % Do we run the pure-state simplification of the hierarchy? Default: false.
options.projective = false; % Do we run the pure-state simplification of the hierarchy? Default: false.
options.forceCompletenessConstraints = true; % Even with projective measurements, this provides tighter bound for Inequality 4
% We use the variant of the hierarchy built around \sqrt{E_{b|y}} operators
options.sqrtRhoHierarchy = false;
options.sqrtEHierarchy = true;

% We can specify yalmip options (e.g., solver, etc.). Otherwise, default solver will be used.
% options.yalmipOptions = sdpsettings('verbose',1,'solver','mosek');
% To run with the levels here (at least for some inequalities) requirs *loads* of memory for mosek, but 
% easily doable with scs (https://github.com/bodono/scs-matlab)
options.yalmipOptions = sdpsettings('verbose',1,'dualize',1,'solver','scs-direct','scs.eps',1e-5,'scs.max_iters',30000);
options.symmetrise = 1; 

% We do this just so code doesns't crash if we comment out some inequalities below and solve specific ones...
score1 = 0; score2 = 0; score3 = 0; score4 = 0; score5 = 0; score6 = 0; score7 = 0;

options.symmetryGenerators = symmetryGenerators1;
[score1, yalmipOut1, Gamma1, Lambda1, Upsilon1] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w1, levels, options);
options.symmetryGenerators = symmetryGenerators2;
[score2, yalmipOut2, Gamma2, Lambda2, Upsilon2] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w2, levels, options);
options.symmetryGenerators = symmetryGenerators3;
[score3, yalmipOut3, Gamma3, Lambda3, Upsilon3] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w3, levels, options);
options.symmetryGenerators = symmetryGenerators4;
[score4, yalmipOut4, Gamma4, Lambda4, Upsilon4] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w4, levels, options);
options.symmetryGenerators = symmetryGenerators5;
[score5, yalmipOut5, Gamma5, Lambda5, Upsilon5] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w5, levels, options);
options.symmetryGenerators = symmetryGenerators6;
[score6, yalmipOut6, Gamma6, Lambda6, Upsilon6] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w6, levels, options);
options.symmetryGenerators = symmetryGenerators7;
[score7, yalmipOut7, Gamma7, Lambda7, Upsilon7] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w7, levels, options);

[NC_bounds; Q_bounds; score1 score2 score3 score4 score5 score6 score7]
