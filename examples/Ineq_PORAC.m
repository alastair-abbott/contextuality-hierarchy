% Inequality from A. Ambainis, M. Banik, A. Chaturvedi, D. Kravchenko and A. Rai, Quantum Inf Process 18, 111 (2019).
% Quantum bounds: d=2: 0.85355, d=3: 0.7777777

% This code should work for arbitrary d (tested for d=2 and d=3)
d = 3;
	
% Specify the number of preparations, measurements and outcomes per measurements
nX = d^2; 
nY = 2; 
nB = d;

% Specify the function to maximise: sum_{b,x,y} w(b,x,y)*p(b|x,y)
w = zeros(nB,nX,nY);
for x1 = 0:d-1
    for x2 = 0:d-1
        for y = 0:1
			x = [x1 x2];
			w(x(y+1)+1, d*x1+x2+1, y+1) = 1/(nX*nY);
        end
    end
end

% Constraints on probability distribution: For PORACs we don't have any.
p_constraints = [];

% Specify the operational equivalences:
% Each operational equivalence is a set (cell) of 2x|S_k| matrices
% The first row specifies the sets (S_k), the second the weights for the corresponding preparations (the \xi_k)
opequiv = cell(d,1);
% For PORACs we have one Preperation OE, which is a partition into d sets of d preparations
for s = 0:d-1
	opequiv{s+1} = zeros(2,d);
	for x1 = 0:d-1
		x2 = mod(s-x1,d);
		opequiv{s+1}(1,x1+1) = fromSeveralBases([x1,x2],[d,d])+1; % write x1x2 in base 10
		opequiv{s+1}(2,x1+1) = 1/d;
	end
end
E_P = {opequiv};
E_M = {}; % No measurement operational equivalences

% Now we specify the levels of the hierarchy we want to use.
% These levels are sufficient to find tight quantum bounds for d=2,3.
% 1: states; 2: measurements, 3: sigma (here we have only one (preparation) operational equivalence)
% E.g., [1 2 2] corresponds to including all rho*M*M, terms in the monomial list \mathcal{S}

% Main moment matrix:
% ====== Level 1 (all) ======
levelsS = {1, 2, 3};
% ====== Level 2 terms ======
levelsS{end+1} = [1, 1];
levelsS{end+1} = [1, 2];
levelsS{end+1} = [1, 3];
levelsS{end+1} = [2, 2];
levelsS{end+1} = [2, 3];
% ====== Level 3 terms ======
levelsS{end+1} = [1, 2, 2];
levelsS{end+1} = [2, 2, 3];

% Localising matrices (L for the operational equivalences, O for the operator positivity):
% ====== Level 1 terms ======
levelsL = {1, 2};
% ====== Level 2 terms ======
levelsL{end+1} = [2, 2];

% Use same monomial list for the localising matrices for operator positivity
levelsO = levelsL;

levels = {levelsS, levelsL, levelsO};

% Specify the symmetries in both the states and the measurements (to be understood as simultaneous with those for the states)
% d^2 preparations rho_0,...,rho_{d^2-1} ~= [1,...,d^2]
% nY*nB measurements M_{0|1},...,M_{d-1|1},M_{0|2},...,M_{d-1|2} ~= [1,...,2d]
stateSymGens = cell(2,1);
measurementSymGens = cell(2,1);
% x1 -> x1 + 1 mod d
stateSymGens{1} = [d+1:d^2 1:d]; 
measurementSymGens{1} = [2:d 1 d+1:2*d]; 
% swap x1 and x2
swapX1X2 = @(x) fromSeveralBases(flip(toSeveralBases(x-1,[d,d])),[d,d])+1;
stateSymGens{2} = arrayfun(swapX1X2,1:d^2);
measurementSymGens{2} = [d+1:2*d 1:d];

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
options.yalmipOptions = sdpsettings('verbose',1,'solver','mosek');

[score, yalmipOut, Gamma, Lambda, Upsilon] = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w, levels, options);