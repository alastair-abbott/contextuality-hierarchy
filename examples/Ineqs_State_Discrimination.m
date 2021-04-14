% Family of noncontextuality inequalites based on state discrimination, from D. Schmid and R. Spekkens, Phys. Rev. X 8, 011015 (2018).
% The inequality is given the Eq. (36) therein, but note, there is a typo in the published version
% (the (c-1) under the square root should be (1-c), as is easily verified).
% The inequalities are paramterised by the error, eps, and "confusability", c, with eps ≤ c ≤ 1 - eps.

% The analytic bound from the paper
analyticBound = @(cc,ee) 1/2*(1+sqrt(1-ee+2*sqrt(ee*(1-ee)*cc*(1-cc))+cc*(2*ee-1)));

stepSize = 0.01;

sHierarchy = [];
sAnalytic = [];
i = 1;
for eps = 0:stepSize:0.5
	for c = eps:stepSize:1-eps-1e-6
		hierarchyVals(i,:) = [c, eps, maximise_SD_inequality(c,eps)];
		analyticVals(i,:) = [c, eps, analyticBound(c,eps)];
		i = i+1;
	end
end

diffs = abs(hierarchyVals(:,3) - analyticVals(:,3));

disp(['Max deviation from analytic bound: ', num2str(max(diffs)), ', average deviation: ', num2str(mean(diffs))]);

function s = maximise_SD_inequality(c,eps)
% For a given eps and c we run the hierarchy to find an upper bound on s

	assert(eps <= c && c <= 1-eps, 'Error: must have eps <= c <= 1-eps');
	
	% Specify the number of preparations, measurements and outcomes per measurements
	% With respect to the notation of the paper:
	% States: 1 -> P_\phi, 2 -> P_\psi, 3 -> P_\bar{\phi}, 4 -> P_\bar{\psi}
	% Measurements: 1 ~ E_{1|1} -> [\phi|M_\phi], 2 ~ E_{2|1} -> [\bar{\phi}|M_\phi]
	%				3 ~ E_{1|2} -> [\psi|M_\psi], 4 ~ E_{2|2} -> [\bar{\psi}|M_\psi]
	%				5 ~ E_{1|3} -> [g_\phi|M_d],  6 ~ E_{2|3} -> [g_\psi|M_d]
	nX = 4; 
	nY = 3;	
	nB = 2;
	dims = [nB, nX, nY];
	
	% Constraints on probability distribution p(b|x,y)
	% The noncontextuality inequality has the following symmetrised structure (cf. Eqs. (30)-(32))
	% p(1|1,2) = p(1|2,1) = p(2|3,2) = p(2|4,1) = c
	% p(1|1,1) = p(1|2,2) = p(2|3,1) = p(2|4,2) = 1 - eps
	% p(1|1,3) = p(2|2,3) = p(2|3,3) = p(1|4,3) (= s)
	%
	% Generically: we give k constraints by specifying a k-by-(nB*nX*nY) matrix A, and a k-element column vector z
	%			   such that A.p = z, where p=[p(1|1,1), p(2|1,1), ..., p(1|2,1), ..., p(nB|nX,nY)]'
	A = [pvecsel(1,1,2,dims); pvecsel(1,2,1,dims); pvecsel(2,3,2,dims); pvecsel(2,4,1,dims)];
	A = [A; pvecsel(1,1,1,dims); pvecsel(1,2,2,dims); pvecsel(2,3,1,dims); pvecsel(2,4,2,dims)];
	A = [A; pvecsel(1,1,3,dims) - pvecsel(2,2,3,dims)];
	A = [A; pvecsel(1,1,3,dims) - pvecsel(2,3,3,dims)];
	A = [A; pvecsel(1,1,3,dims) - pvecsel(1,4,3,dims)];
	
	z = [c; c; c; c; 1-eps; 1-eps; 1-eps; 1-eps; 0; 0; 0];
	p_constraints = {A,z};
	
	% Specifies the score function for the task, so that score = sum_{b,x,y} w(b,x,y)*p(b|x,y)
	% Here we just want to maximise s = p(1|1,3), but we symmetrise it so it works better with the symmetries we impose below
	w = zeros(nB,nX,nY);
	w(1,1,3) = 1/2;
	w(2,2,3) = 1/2;

	% Specify the operational equivalences:
	% Each operational equivalence is a set (cell) of 2x|S_k| matrices
	% The first row specifies the sets (S_k), the second the weights for the corresponding preparations (the \xi_k)
	opequiv = cell(2,1);
	opequiv{1} = [1 3; 1/2 1/2];
	opequiv{2} = [2 4; 1/2 1/2];

	E_P = {opequiv};
	E_M = {}; % No measurement operational equivalences

	% Now we specify the levels of the hierarchy we want to use
	% 1: states; 2: measurements, 3: sigma (here we have only one (preparation) operational equivalence)
	% E.g., [1 2 2] corresponds to including all rho*M*M, terms in the monomial list \mathcal{S}

	% Main moment matrix:
	% ====== Level 1 ======
	levelsS = {1, 2, 3};
	% ====== Level 2 ======
	levelsS{end+1} = [1, 2];
	levelsS{end+1} = [2, 3];

	% Localising matrices (L for the operational equivalences, O for the operator positivity):
	% ====== Level 1 ======
	levelsL = {2};

	% Use same monomial list for the localising matrices for operator positivity
	levelsO = levelsL;

	levels = {levelsS, levelsL, levelsO};
	
	% Specify the symmetries in both the states and the measurements (to be understood as simultaneous with those for the states)
	stateSymGens = cell(1,1);
	measurementSymGens = cell(1,1);
	% Swap \psi <-> \phi, \bar{\psi} <-> \bar{\phi}, M_{b|\psi} <-> M_{b|\phi}, M_{g_\psi|d} <-> M_{g_\phi|d}
	stateSymGens{1} = [2 1 4 3];
	measurementSymGens{1} = [3 4 1 2 6 5];

	symmetryGenerators = {stateSymGens, measurementSymGens};

	% Options for hierarchy
	options = struct();
	options.symmetrise = 0; % The problem here is small enough that no real advantage from symmetrising, but turn can enable for fun ;-)
	options.symmetryGenerators = symmetryGenerators;

	% options.classical = true; % Do we want the classical (commuting) relaxation? Default: false.
	% options.pure = true; % Do we run the pure-state simplification of the hierarchy? Default: false.

	% We can specify yalmip options (e.g., solver, etc.). Otherwise, default solver will be used.
	% options.yalmipOptions = sdpsettings('verbose',1,'solver','mosek');

	s = boundContextualCorrelations([nX, nY, nB], E_P, E_M, p_constraints, w, levels, options);

end

function v = pvecsel(b,x,y,dims)
% Returns an array v such that v.p = p(b|x,y), where p=[p(1|1,1), p(2|1,1), ..., p(1|2,1), ..., p(nB|nX,nY)]
% In other words, v is 0 everyone except at the element corresponding to p(b|x,y), where it is 1
% We expect dims = [nB, nX, nY]
	v = zeros(dims);
	v(b,x,y) = 1;
	v = v(:).';
end