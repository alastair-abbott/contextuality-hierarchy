% Noncontextuality inequality from parity oblivious multiplexity, 
% from R. W. Spekkens, D. H. Buzacott, A. J. Keehn, Ben Toner, G. J. Pryde, Phys. Rev. Lett. 102, 010401 (2009).
% Here we implement the n=2 and n=3 cases.

%% Parity oblivious multiplexing: n = 2
% Noncontextual bound: 3/4
% Optimal quantum violation: cos^2(pi/8) = 0.85355

% Specify the number of preparations, measurements and outcomes per measurements
nX = 4; 
nY = 2; 
nB = 2;

% Specifies the score function for the task, so that score = sum_{b,x,y} w(b,x,y)*p(b|x,y)
w = zeros(nB,nX,nY);
for x1 = 0:1
    for x2 = 0:1
        for y = 0:1
            x = [x1 x2];
			w(x(y+1)+1, x1*2 + x2 + 1, y+1) = 1/(nX*nY);
        end
    end
end

% Constraints on probability distribution: For this inequality we don't have any.
p_constraints = [];

% Specify the operational equivalences:
% Each operational equivalence is a set (cell) of 2x|S_k| matrices
% The first row specifies the sets (S_k), the second the weights for the corresponding preparations (the \xi_k)
opequiv = cell(2,1);
opequiv{1} = [1 4; 1/2 1/2];
opequiv{2} = [2 3; 1/2 1/2];

E_P = {opequiv};
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
levelsS{end+1} = [2, 1];
levelsS{end+1} = [2, 2];
levelsS{end+1} = [2, 3];
levelsS{end+1} = [3, 3];

% Localising matrices (L for the operational equivalences, O for the operator positivity):
% ====== Level 1 ======
levelsL = {1, 2};

% If we want to use the classical restriction to re-obtain the noncontextual bound, we need some extra levels
% levelsS{end+1} = [1, 2, 2];
% levelsS{end+1} = [2, 2, 3];
% levelsL{end+1} = [2, 2];

% Use same monomial list for the localising matrices for operator positivity
levelsO = levelsL;

levels = {levelsS, levelsL, levelsO};

% Specify the symmetries in both the states and the measurements (to be understood as simultaneous with those for the states)
% nX preparations rho_0,...,rho_{nX-1} ~= [1,...,nX]
% nY*nB measurements M_{0|0},...,M_{nB-1|nY-1} ~= [1,...,nY*nB]
stateSymGens = cell(2,1);
measurementSymGens = cell(2,1);
% x0 -> x0 + 1 mod 2 (where x=x0x1 in base 2)
stateSymGens{1} = [3 4 1 2];
measurementSymGens{1} = [2 1 3 4]; % nY*nB measurements; we handle the offset internally, here we just label measurements 1:6
% swap x0 and x1
stateSymGens{2} = [1 3 2 4];
measurementSymGens{2} = [3 4 1 2];

symmetryGenerators = {stateSymGens, measurementSymGens};

% Options for hierarchy
options = struct();
options.verbose = true;
% Symmetrisation: 
%	0: no symmetrisation (deault); no need to specify symmetryGenerators in that case
%	1: invariant moment matrix; need to specify symmetryGenerators
%	2: block-diagonalisation (experimental - not recommended); requires symmetryGenerators. Smaller SDP but slower pre-processing
% options.symmetrise = 0; % Here we can easily solve without symmetrising, so no need to specify generators
% options.symmetryGenerators = symmetryGenerators;

% options.classical = true; % Do we want the classical (commuting) relaxation? Default: false.
% options.pure = true; % Do we run the pure-state simplification of the hierarchy? Default: false.

% We can specify yalmip options (e.g., solver, etc.). Otherwise, default solver will be used.
% options.yalmipOptions = sdpsettings('verbose',1,'solver','mosek');

[score_2, yalmipOut_2, Gamma_2, Lambda_2, Upsilon_2] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w, levels, options);

%% Parity oblivious multiplexing: n = 3 
% Noncontextual bound: 2/3
% Optimal quantum violation: 0.788675

% Specify the number of preparations, measurements and outcomes per measurements
nX = 8; 
nY = 3; 
nB = 2;

% Score function: sum_{b,x,y} w(b,x,y)*p(b|x,y)
w = zeros(nB,nX,nY);
for x1 = 0:1
    for x2 = 0:1
		for x3 = 0:1
			for y = 0:2
				x = [x1 x2 x3];
				w(x(y+1)+1, x1*4 + x2*2 + x3 + 1, y+1) = 1/(nX*nY);
			end
		end
    end
end

% Constraints on probability distribution: For this inequality we don't have any.
p_constraints = [];

% Specify the operational equivalences (now we have 4):
E_P = cell(4,1);
E_P{1} = cell(2,1);
E_P{1}{1} = [1 2 7 8; 1/4 1/4 1/4 1/4];
E_P{1}{2} = [3 4 5 6; 1/4 1/4 1/4 1/4];
E_P{2} = cell(2,1);
E_P{2}{1} = [1 3 6 8; 1/4 1/4 1/4 1/4];
E_P{2}{2} = [2 4 5 7; 1/4 1/4 1/4 1/4];
E_P{3} = cell(2,1);
E_P{3}{1} = [1 4 5 8; 1/4 1/4 1/4 1/4];
E_P{3}{2} = [2 3 6 7; 1/4 1/4 1/4 1/4];
E_P{4} = cell(2,1);
E_P{4}{1} = [1 4 6 7; 1/4 1/4 1/4 1/4];
E_P{4}{2} = [2 3 5 8; 1/4 1/4 1/4 1/4];

E_M = {}; % No measurement operational equivalences

% Now we specify the levels of the hierarchy we want to use
% Main moment matrix:
% ====== Level 1 ======
levelsS = {1, 2, 3};
% ====== Level 2 ======
levelsS{end+1} = [1, 1];
levelsS{end+1} = [1, 2];
levelsS{end+1} = [1, 3];
levelsS{end+1} = [2, 1];
levelsS{end+1} = [2, 2];
levelsS{end+1} = [2, 3];
levelsS{end+1} = [3, 3];

% Localising matrices (L for the operational equivalences, O for the operator positivity):
% ====== Level 1 ======
levelsL = {1, 2};

% Use same monomial list for the localising matrices for operator positivity
levelsO = levelsL;

levels = {levelsS, levelsL, levelsO};

% Specify the symmetries
stateSymGens = cell(3,1);
measurementSymGens = cell(3,1);
% x1 -> x1 + 1 mod 2 (where x=x1x2x3 in base 2)
stateSymGens{1} = [5:8 1:4];
measurementSymGens{1} = [2 1 3:6]; % nY*nB measurements; we handle the offset internally, here we just label measurements 1:6
% swap x1 and x2
stateSymGens{2} = [1 2 5 6 3 4 7 8];
measurementSymGens{2} = [3 4 1 2 5 6];
% cycle x1->x2->x3->x1
stateSymGens{3} = [1 5 2 6 3 7 4 8];
measurementSymGens{3} = [3:6 1:2];

symmetryGenerators = {stateSymGens, measurementSymGens};

% Options for hierarchy
options = struct();
options.verbose = true;
% Symmetrisation: 
%	0: no symmetrisation (deault); no need to specify symmetryGenerators in that case
%	1: invariant moment matrix; need to specify symmetryGenerators
%	2: block-diagonalisation (experimental - not recommended); requires symmetryGenerators. Smaller SDP but slower pre-processing
options.symmetrise = 1; % Unlike the n=2 case, symmetrising helps a lot here, although still solvable withouot
options.symmetryGenerators = symmetryGenerators;

% options.pure = true; % Do we run the pure-state simplification of the hierarchy? Default: false.

% We can specify yalmip options (e.g., solver, etc.). Otherwise, default solver will be used.
% options.yalmipOptions = sdpsettings('verbose',1,'solver','mosek');

[score_3, yalmipOut_3, Gamma_3, Lambda_3, Upsilon_3] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w, levels, options);
