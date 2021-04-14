% Simple inequality from A. Hameedi, A. Tavakoli, B. Marques, and M. Bourennane, Phys. Rev. Lett. 119, 220402 (2017)
% c.f. Eq. (9) therein. The inequality is equivalent to the CGLMP Bell inequality, with noncontextual bound of 1/2
% and optimal quantum violation giving 0.7287.

% Specify the number of preparations, measurements and outcomes per measurements
nX = 6; 
nY = 2; 
nB = 3;

% Specifies the score function for the task, so that score = sum_{b,x,y} w(b,x,y)*p(b|x,y)
w = zeros(nB,nX,nY);
for x1 = 0:1
    for x2 = 0:2
        for y = 0:1
			for m = 0:1
				x = x1*3 + x2;
				T = mod(x2 - (-1)^(x1+y+m)*m - x1*y, 3);
				w(T+1, x+1, y+1) = w(T+1, x1*3 + x2 + 1, y+1) + (-1)^m/(nX*nY);
			end
        end
    end
end

% Constraints on probability distribution: For this inequality we don't have any.
p_constraints = [];

% Specify the operational equivalences:
% Each operational equivalence is a set (cell) of 2x|S_k| matrices
% The first row specifies the sets (S_k), the second the weights for the corresponding preparations (the \xi_k)
opequiv = cell(2,1);
opequiv{1} = [1 2 3; 1/3 1/3 1/3];
opequiv{2} = [4 5 6; 1/3 1/3 1/3];

E_P = {opequiv};
E_M = {}; % No measurement operational equivalences

% Now we specify the levels of the hierarchy we want to use
% These levels are sufficient to find the tight quantum bound here
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
levelsS{end+1} = [1, 2, 2];
levelsS{end+1} = [2, 2, 3];

% Localising matrices (L for the operational equivalences, O for the operator positivity):
% ====== Level 1 ======
levelsL = {1, 2};
% ====== Level 2 ======
levelsL{end+1} = [2, 2];

% Use same monomial list for the localising matrices for operator positivity
levelsO = levelsL;

levels = {levelsS, levelsL, levelsO};

% Specify the symmetries in both the states and the measurements (to be understood as simultaneous with those for the states)
% nX preparations rho_0,...,rho_{nX-1} ~= [1,...,nX]
% nY*nB measurements M_{1|1},...,M_{nB|nY} ~= [1,...,nY*nB]
stateSymGens = cell(1,1);
measurementSymGens = cell(1,1);
% x2 -> x2 + 1 mod 3
stateSymGens{1} = [2 3 1 5 6 4];
measurementSymGens{1} = [2 3 1 5 6 4];

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
% options.pure = true; % Do we run the pure-state simplification of the hierarchy? Default: false.

% We can specify yalmip options (e.g., solver, etc.). Otherwise, default solver will be used.
% options.yalmipOptions = sdpsettings('verbose',1,'solver','mosek');

[score, yalmipOut, Gamma, Lambda, Upsilon] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w, levels, options);