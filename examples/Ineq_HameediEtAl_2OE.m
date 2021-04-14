% Three outcome noncontextuality inequality with multiple operation equivalences from the Appendix of 
% A. Hameedi, A. Tavakoli, B. Marques, and M. Bourennane, Phys. Rev. Lett. 119, 220402 (2017).

% Here we take the n=2, d=3 case, for which the noncontextual bound is 2/3.
% Specify the number of preparations, measurements and outcomes per measurements
nX = 9;
nY = 2; 
nB = 3; 

% Specifies the score function for the task, so that score = sum_{b,x,y} w(b,x,y)*p(b|x,y)
w = zeros(nB,nX,nY);
for x1 = 0:2
    for x2 = 0:2
        for y = 0:1
            x = [x1 x2];
			w(x(y+1)+1, x1*3 + x2 + 1, y+1) = 1/(nX*nY);
        end
    end
end

% Constraints on probability distribution: For this inequality we don't have any.
p_constraints = [];

% Specify the operational equivalences:
% Each operational equivalence is a set (cell) of 2x|S_k| matrices
% The first row specifies the sets (S_k), the second the weights for the corresponding preparations (the \xi_k)
E_P = cell(2,1);
E_P{1} = cell(3,1);
E_P{1}{1} = [1 6 8; 1/3 1/3 1/3];
E_P{1}{2} = [2 4 9; 1/3 1/3 1/3];
E_P{1}{3} = [3 5 7; 1/3 1/3 1/3];
E_P{2} = cell(3,1);
E_P{2}{1} = [1 5 9; 1/3 1/3 1/3];
E_P{2}{2} = [3 4 8; 1/3 1/3 1/3];
E_P{2}{3} = [2 6 7; 1/3 1/3 1/3];

E_M = {}; % No measurement operational equivalences

% Now we specify the levels of the hierarchy we want to use
% 1: states; 2: measurements, 3: sigma (here we have only one (preparation) operational equivalence)
% E.g., [1 2 2] corresponds to including all rho*M*M, terms in the monomial list \mathcal{S}

% Main moment matrix:
% ====== Level 1 ======
levelsS = {1, 2, 3};
% ====== Level 2 ======
levelsS{end+1} = [1, 1];
levelsS{end+1} = [1, 2];
levelsS{end+1} = [1, 3];
levelsS{end+1} = [2, 2];
levelsS{end+1} = [2, 3];
% ====== Level 3 ======
levelsS{end+1} = [1, 1, 1];
levelsS{end+1} = [1, 1, 3];
levelsS{end+1} = [1, 2, 1];
levelsS{end+1} = [1, 2, 2];
levelsS{end+1} = [1, 2, 3];
levelsS{end+1} = [2, 2, 3];
levelsS{end+1} = [2, 3, 2];

% Localising matrices (L for the operational equivalences, O for the operator positivity):
% ====== Level 1 ======
levelsL = {1, 2};
% ====== Level 2 ======
levelsL{end+1} = [1, 1];
levelsL{end+1} = [2, 1];
levelsL{end+1} = [2, 2];

% % If we want a better upper bound for quantum models (i.e., not with pure states), 
% % we need more levels (and set the pure option to false below)
%
% % Main moment matrix:
% % ====== Level 1 ======
% levelsS = {1, 2, 3};
% % ====== Level 2 ======
% levelsS{end+1} = [1, 1];
% levelsS{end+1} = [1, 2];
% levelsS{end+1} = [1, 3];
% levelsS{end+1} = [2, 1];
% levelsS{end+1} = [2, 2];
% levelsS{end+1} = [2, 3];
% levelsS{end+1} = [3, 1];
% levelsS{end+1} = [3, 2];
% levelsS{end+1} = [3, 3];
% % ====== Level 3 ======
% levelsS{end+1} = [1, 2, 2];
% levelsS{end+1} = [2, 2, 2];
% levelsS{end+1} = [2, 2, 3];
% % ====== Level 4 ======
% levelsS{end+1} = [1, 2, 2, 2];
% levelsS{end+1} = [2, 2, 2, 3];
% 
% % The operational equivalence localising matrices
% % ====== Level 1 ======
% levelsL = {1, 2};
% % ====== Level 2 ======
% levelsL{end+1} = [2, 2];
% % ====== Level 2 ======
% levelsL{end+1} = [2, 2, 2];

% Use same monomial list for the localising matrices for operator positivity
levelsO = levelsL;

levels = {levelsS, levelsL, levelsO};

% Specify the symmetries in both the states and the measurements (to be understood as simultaneous with those for the states)
% nX preparations rho_0,...,rho_{nX-1} ~= [1,...,nX]
% nY*nB measurements M_{1|1},...,M_{nB|nY} ~= [1,...,nY*nB]
stateSymGens = cell(2,1);
measurementSymGens = cell(2,1);
% x0 -> x0 + 1 mod 3 (where x=x0x1 in base 3)
stateSymGens{1} = [4:9 1:3];
measurementSymGens{1} = [2 3 1 4:6];
% swap x0 and x1
stateSymGens{2} = [1 4 7 2 5 8 3 6 9];
measurementSymGens{2} = [4:6 1:3];

symmetryGenerators = {stateSymGens, measurementSymGens};

% Options for hierarchy
options = struct();
options.verbose = true;
% Symmetrisation: 
%	0: no symmetrisation (deault); no need to specify symmetryGenerators in that case
%	1: invariant moment matrix; need to specify symmetryGenerators
%	2: block-diagonalisation (experimental - not recommended); requires symmetryGenerators. Smaller SDP but slower pre-processing
options.symmetrise = 1; 
options.symmetryGenerators = symmetryGenerators;

% options.classical = true; % Do we want the classical (commuting) relaxation? Default: false.
options.pure = true; % Do we run the pure-state simplification of the hierarchy? Default: false.

% We can specify yalmip options (e.g., solver, etc.). Otherwise, default solver will be used.
options.yalmipOptions = sdpsettings('verbose',1,'solver','mosek');
% options.yalmipOptions = sdpsettings('verbose',1,'dualize',1,'solver','scs-direct','scs.eps',1e-6,'scs.max_iters',50000);

[score, yalmipOut, Gamma, Lambda, Upsilon] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w, levels, options);