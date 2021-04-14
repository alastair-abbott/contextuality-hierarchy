% Noncontextuality inequality from parity oblivious multiplexity, 
% from R. W. Spekkens, D. H. Buzacott, A. J. Keehn, Ben Toner, G. J. Pryde, Phys. Rev. Lett. 102, 010401 (2009).
% Here we bound the simulation cost for n = 2

%% Specify the number of preparations, measurements and outcomes per measurements
nX = 4; 
nY = 2; 
nB = 2;

% Specifies the score function for the task, so that score = sum_{b,x,y} w(b,x,y)*p(b|x,y)
% We want to calculate the cost of simulating a specific score
w = zeros(nB,nX,nY);
for x1 = 0:1
    for x2 = 0:1
        for y = 0:1
            x = [x1 x2];
			w(x(y+1)+1, x1*2 + x2 + 1, y+1) = 1/(nX*nY);
        end
    end
end
w = w(:);

% Specify the operational equivalences:
% Each operational equivalence is a set (cell) of 2x|S_k| matrices
% The first row specifies the sets (S_k), the second the weights for the corresponding preparations (the \xi_k)
opequiv = cell(2,1);
opequiv{1} = [1 4; 1/2 1/2];
opequiv{2} = [2 3; 1/2 1/2];

E_P = {opequiv};

% Now we specify the levels of the hierarchy we want to use
% These levels are sufficient to obtain the tight quantum simulation cost
% The classical simulation cost can be recovered without the Level 4 terms in S or level 3 terms in L and O.

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
% ====== Level 3 ======
levelsS{end+1} = [1, 2, 2];
levelsS{end+1} = [2, 2, 2];
levelsS{end+1} = [2, 2, 3];
% ====== Level 4 ======
levelsS{end+1} = [1, 2, 2, 2];
levelsS{end+1} = [2, 2, 2, 3];

% Localising matrices (L for the operational equivalences, O for the operator positivity):
% ====== Level 1 ======
levelsL = {1, 2};
% ====== Level 2 ======
levelsL{end+1} = [2, 2];
% ====== Level 3 ======
levelsL{end+1} = [2, 2, 2];

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
options.verbose = false;
% Symmetrisation: 
%	0: no symmetrisation (deault); no need to specify symmetryGenerators in that case
%	1: invariant moment matrix; need to specify symmetryGenerators
%	2: block-diagonalisation (experimental - not recommended); requires symmetryGenerators. Smaller SDP but slower pre-processing
options.symmetrise = 1; % Here we can easily solve without symmetrising, so no need to specify generators
options.symmetryGenerators = symmetryGenerators;

options.classical = false; % Bound the classical simulation cost instead? Default: false
% options.pure = true; % Do we run the pure-state simplification of the hierarchy? Default: false.

% We can specify yalmip options (e.g., solver, etc.). Otherwise, default solver will be used.
% options.yalmipOptions = sdpsettings('verbose',1,'solver','mosek');

numSamples = 26;
% score = 1/2*(1+1/sqrt(2)):1/2*(1-1/sqrt(2))/(numSamples-1):1;
score = 0.75:0.25/(numSamples-1):1;
simulationCost = zeros(1,numSamples);
G = zeros(1,numSamples);

for i = 1:numSamples
	display(['Starting sample number ', num2str(i)]);
	% We calculate the simulation cost under the constraint that w.'*p = dot(w,p) = score(i)
	p_constraints = {w.', score(i)};
	[simulationCost(i), G(i)] = boundSimulationCost([nX, nY, nB], E_P, p_constraints, levels, options);
end