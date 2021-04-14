function [cost, G_out, yalmipOut, GammaOpt, LambdaOpt, UpsilonOpt] = boundSimulationCost(n_xyb, E_P, p_constraints, levels, options)
%boundSimulationCost Bounds the simulation cost of preparation contextual correlations satisfying given constraints using the SDP hierarchy
%                    described in the paper: "A. Tavakoli, E. Zambrini Cruzeiro, R. Uola, A. A. Abbott, 
%                    "Bounding and simulating contextual correlations in quantum theory", arXiv:2010.04751 [quant-ph] (2020)."
%   Output:
%       - cost: the simulation cost in bits
%       - G_out: the Guessing probability, so that cost = log2(K) + log2(G_out)
%       - yalmipOut: the yalmip output, containing notably the solver messages and dimacs error metrics
%       - GammaOpt, LambdaOpt, UpsilonOpt: the matrices corresponding to the optimal SDP solution
%   Arguments:
%		- n_xyb = [nX, nY, nB]: the number of states, measurements and outcomes, respectively.
%                 * For the rest of the options, the states are assumed to be labelled 1,...,nX
%                  and measurement effects 1,...,nB*nY, with the ordering E_{1|1},...,E_{nB|1},E_{1|2},...,E_{nB|nY}
%		- E_P: A cell array such that E_P{r} is a 1xK_r cell-array specifying the r-th preparation operational equivalence,
%                 i.e., E_P{r}{k} is a 2 row matrix; the first row labelling the states x in the set S_k, 
%                 the second row the corresponding probabilities \xi_k(x)
%		- p_constraints = {A,z}: constraints on the distribution to be simulated, so that the probability distribution p satisfies A.p = z;
%				  can be empty ({} or []) if no constraints to impose (but then simulation cost is trivial).
%                 p is taken to be the probability vector (p(1|1,1),...,p(nB|1,1),p(1|2,1),...).
%		- levels = {levelsS, levelsL, levelsO}: 3 cell arrays specifying the hierarchy levels for the moment matrix and respective localising matrices.
%                 each level is an array over {1,2,3}, where 1 specifies states, 2 measurements, and 3 sigma/tau operators (used to enforce operational equivalences)
%       - options: an optional structure specifying further options for the sdp hierarchy:
%                 > verbose (default: false): Print out helpful information while running hierarchy.
%                 > classical (default: false): Treat all operators as commuting.
%                 > pure (default: false): Treat all states as being pure.
%                 > symmetrise (default: 0): Specify whether the moment matrix should be symmetrised.
%                       - 0: No symmetrisation applied.
%                       - 1: Make the moment matrix invariant under symmetries. Requires replab and specifying symmetryGenerators (see below).
%                       - 2: Block diagonalise the moment matrix. Requires replab and specifying symmetryGenerators.
%                 > symmetryGenerators = {stateSymGens, measurementSymGens} (required if symmetrise = 1,2): 
%                       * stateSymGens: a cell array of permutations on the states {1,...,nX} generating the symmetries
%                       * measurementSymGens: a cell array of the corresponding permutations on the measurements {1,...,nY*nB}
%                 > yalmipOptions (defaults to default solver and verbose=0). Options to pass to yalmip to specify, e.g., the solver or verbose output
%                 > forceCompletenessConstraints (default: false): Impose completeness relation even if measurements projective

	%% Check and process the input
	t_total = tic;

	assert(length(n_xyb) == 3 && all(n_xyb > 0), 'Error: n_xyb should be a vector [nX, nY, nB] with nX,nY,nB > 0.');
	nX = n_xyb(1);
	nY = n_xyb(2);
	nB = n_xyb(3);
	
	assert(isempty(E_P) || iscell(E_P), 'Error: operational equivalences not provided in required form.');
	R = length(E_P);
	assert(R == 1, 'Error: simulation cost only defined when a single operational equivalence present (for now).');
	
	% Check the validty of the OEs
	for r = 1:R
		opequiv = E_P{r};
		assert(iscell(opequiv), 'Error: preparation operational equivalence not provided in required form.');
		K = length(opequiv);
		statesInOpequiv = [];
		for k = 1:K
			statesInOpequiv = [statesInOpequiv, opequiv{k}(1,:)];
			assert(all(opequiv{k}(2,:) > 0) && sum(opequiv{k}(2,:)) == 1, ['Error: Weights for preparation operational equivalence ', num2str(r), ', set ', num2str(k), ' not a valid probability distribution.']);
		end
		assert(isequal(sort(statesInOpequiv),1:nX), ['Error: Preparation operational equivalence ', num2str(r), ' not given over a valid partition of [nX].']);
	end
	
	% Check constraints specified properly
	if ~isempty(p_constraints)
		assert(iscell(p_constraints) && length(p_constraints) == 2, 'Error: Constraints should be specified as a pair {A,z}.');
		
		A = p_constraints{1};
		z = p_constraints{2};
		
		assert(nB*nX*nY == size(A,2), 'Error: Constraint matrix A of incorrect dimensions.');
		assert(length(z) == size(A,1), 'Error: Constraint matrix A and vector z should have same number of rows.');
	end
		
	assert(iscell(levels) && isequal(size(levels),[1 3]), 'Error: levels should be a cell {levelsS, levelsL, levelsO} of the levels for each list of moments.');
	levelsS = levels{1};
	levelsL = levels{2};
	levelsO = levels{3};
	
	% Check levels specified correctly
	assert(iscell(levelsS) && size(levelsS,1) == 1, 'Error: levelsS should be a cell indicating which levels of the hierarchy to use.');
	assert(iscell(levelsL) && size(levelsL,1) == 1, 'Error: levelsL should be a cell indicating which levels of the hierarchy to use.');
	assert(iscell(levelsO) && size(levelsO,1) == 1, 'Error: levelsO should be a cell indicating which levels of the hierarchy to use.');
	for i = 1:length(levelsS)
		assert(all(ismember(unique(levelsS{i}),[1 2 3])), 'Error: Each level should be specified of a list composed of 1 (states), 2 (measurements), 3 (operational equivalence operators, i.e. sigma/tau).')
	end
	for i = 1:length(levelsL)
		assert(all(ismember(unique(levelsL{i}),[1 2 3])), 'Error: Each level should be specified of a list composed of 1 (states), 2 (measurements), 3 (operational equivalence operators, i.e. sigma/tau).')
	end
	for i = 1:length(levelsO)
		assert(all(ismember(unique(levelsO{i}),[1 2 3])), 'Error: Each level should be specified of a list composed of 1 (states), 2 (measurements), 3 (operational equivalence operators, i.e. sigma/tau).')
	end
	
	% Default options
	dfoptions.verbose = false;
	dfoptions.symmetrise = 0;
	dfoptions.classical = false;
	dfoptions.pure = false;
	dfoptions.forceCompletenessConstraints = false;
	dfoptions.yalmipOptions = sdpsettings('verbose',0);
	
	if ~exist('options','var')
		options = dfoptions;
	else
		assert(isstruct(options),'Error: options should be a structure.');
		% Set options to default values
		optionNames = fieldnames(dfoptions);
		for i = 1:length(optionNames)
			if ~isfield(options,optionNames{i})
				options.(optionNames{i}) = dfoptions.(optionNames{i});
			end
		end
		% Check options are valid
		assert(ismember(options.symmetrise,[0,1,2]),'Error: invalid symmetry option.');
		if ismember(options.symmetrise,[1,2])
			assert(isfield(options,'symmetryGenerators'),'Error: no symmetry generators provided.')
			% Check symmetries are given as valid permutations
			symS = options.symmetryGenerators{1};
			symM = options.symmetryGenerators{2};
			assert(length(symS) == length(symM),'Error: Every state symmetry must have a corresponding measurement symmetry specified.');
			warning_nonInvariantConstr = false;
			for i = 1:length(symS)
				assert(isequal(sort(symS{i}),1:nX) && isequal(sort(symM{i}),1:nY*nB),['Error: Symmetry number ', num2str(i), ' not specified as a valid permutation of states and measurements.'])
				% Also check that the constraints are indeed invariant under the symmetries
				pindices = reshape(1:nB*nX*nY,[nB,nX,nY]);
				psym = zeros(nB,nX,nY);
				for x = 1:nX
					for y = 1:nY
						for b = 1:nB
							yb = fromSeveralBases([y-1,b-1], [nY,nB]) + 1;
							newYB = toSeveralBases(symM{i}(yb)-1,[nY,nB]) + [1,1];
							
							psym(newYB(2),symS{i}(x),newYB(1)) = pindices(b,x,y);
						end
					end
				end
				
				if ~isempty(p_constraints)
					psym = psym(:);
					Asym = A(:,psym);
					% Check if A and Asym describe same constraints
					if max(abs(A\z - Asym\z)) > 1e-9
						warning_nonInvariantConstr = true;
					end
				end
			end
			if warning_nonInvariantConstr
				disp('!!! Warning: constraints may not be invariant under supplied symmetries: please double check (only rudimentary check performed) !!!');
			end
		end
		assert(isstruct(options.yalmipOptions),'Error: invalid yalmip options - should be a sturcture.');
	end
	
	symmetrise = options.symmetrise;
	classical = options.classical;
	pure = options.pure;
	projective = true; % Here we are simulating preparation noncontextuality, so can always assume projectivity
	forceCompletenessConstraints = options.forceCompletenessConstraints;
	yalmipOptions = options.yalmipOptions;
	
	% Symmetrisation requires replab to be installed and in path
	if symmetrise > 0
		try
			replab_init;
		catch
			disp('Error initialising replab. Check replab installed correctly and present in path.');
			disp('Continuing without symmetrisation...');
			symmetrise = 0;
		end
		if (symmetrise == 1 && ~exist('+replab/+cvar/symmetrizeIndexMatrix')) || ...
				(symmetrise == 2 && ~exist('+replab/CommutantVar'))
			disp('Replab version does not support required symmetrisation approach. Please update replab to a recent develop build (tested with commit 6f2fdd).');
			disp('Continuing without symmetrisation...');
			symmetrise = 0;
		end
	end
	
	if R == 0 && ~pure
		vdisp('No preparation operational equivalences: assuming pure states');
		pure = true;
	end
	
	%% Set up some generic initial stuff

	% The indices for each operator we're using. Index 0 is reserved for the identity.
	states = 1:nX;
	measurements = nX+1:nX+nY*nB;
	sigmas = measurements(end)+1:measurements(end)+R;
	numOpTypes = [nX, nY*nB, R];
	numOps = sum(numOpTypes);
	opOffsets = [0, nX, nX+nY*nB]; % offsets needed to access corresponding operator type

	isMeasurement = @(index) measurements(1) <= index && index <= measurements(end);

	% Function to determine {b,y} given a measurement, and vice versa
	measurementLabels = @(measurement) toSeveralBases(measurement-measurements(1),[nY,nB])+[1,1];
	measurementFromLabels = @(labels) fromSeveralBases(labels-[1,1], [nY,nB]) + measurements(1);

	%% Construct the list of monomials for the main moment matrix
	% Each monomial is a list [s_1, ..., s_maxLevelS] of indices specifying the corresponding operators
	% The list is right-padded with the identity so that all monomials have same order
	% E.g., [3 4 2 0] corresponding to s_3*s_4*s_2 (with the final identity (0) omitted

	maxLevelS = max(cellfun('length',levelsS));

	% Create the monomial list
	S = zeros(1, maxLevelS); % First we need the identity

	% now loop through each family of monomials we're using
	for i = 1:length(levelsS)
		S0 = toSeveralBases([0:prod(numOpTypes(levelsS{i}))-1]', numOpTypes(levelsS{i})) + opOffsets(levelsS{i}) + 1;
		% now pad to the left with the identity and add to S
		S0 = [zeros(size(S0,1),maxLevelS-size(S0,2)), S0];
		S = [S; S0];
	end
	clear S0;

	numMonomialsS = size(S,1);
	
	% We also make a lookup table to find index of a given monomial in S (much faster than actually checking)
	% We interpret each monomial as an integer in base numOps+1 then use this integer to look in a sparse matrix for the index
	monomialIndexSLookup = sparse((numOps+1)^(maxLevelS),1);
	for i = 1:numMonomialsS
		monomialID = fromSeveralBases(S(i,:),(numOps+1)*ones(1,maxLevelS)) + 1;
		monomialIndexSLookup(monomialID) = i;
	end
	% Define a lookup function that checks if a table of monomials (one per row) are contained in GammaMonomials.
	% Returns the (linear) indices.
	monomialIndexInS = @(monomials) monomialIndices(monomials, monomialIndexSLookup, numOps);
	
	vdisp(['Moment matrix size: ', num2str(numMonomialsS)]);

	%% Construct the list of monomials for the localising matrices for the operational equivalences

	maxLevelL = max(cellfun('length',levelsL));
	assert(maxLevelL < maxLevelS,'Error: Main moment matrix level not large enough to include all operational equivalence localising matrix terms.');

	% Create the monomial list
	L = zeros(1, maxLevelL); % First we need the identity

	% now loop through each family of monomials we're using
	for i = 1:length(levelsL)
		L0 = toSeveralBases([0:prod(numOpTypes(levelsL{i}))-1]', numOpTypes(levelsL{i})) + opOffsets(levelsL{i}) + 1;
		% now pad to the left with the identity and add to L
		L0 = [zeros(size(L0,1),maxLevelL-size(L0,2)), L0];
		L = [L; L0];
	end
	clear L0;

	numMonomialsL = size(L,1);
	vdisp(['Localising matrix size for operational equivalences: ', num2str(numMonomialsL)]);

	%% Construct the list of monomials for the localising matrix for ensuring the states and measurements are PSD

	if ~pure || ~projective % Only need this if we're not doing pure & projective hierarchy
		maxLevelO = max(cellfun('length',levelsO));
		assert(maxLevelO < maxLevelS,'Error: Main moment matrix level not large enough to include all mixed state localising matrix terms.');

		% Create the monomial list
		O = zeros(1, maxLevelO); % First we need the identity

		% now loop through each family of monomials we're using
		for i = 1:length(levelsO)
			O0 = toSeveralBases([0:prod(numOpTypes(levelsO{i}))-1]', numOpTypes(levelsO{i})) + opOffsets(levelsO{i}) + 1;
			% now pad to the left with the identity and add to O
			O0 = [zeros(size(O0,1),maxLevelO-size(O0,2)), O0];
			O = [O; O0];
		end
		clear O0;

		numMonomialsO = size(O,1);
		vdisp(['Localising matrix size for state/measurement positivity: ', num2str(numMonomialsO)]);
	end
	
	%% Generate and simplify the monomial lists for each element of the moment matrix

	t = tic;

	GammaMonomials = zeros(2*maxLevelS,numMonomialsS,numMonomialsS);
	
	% We also make a lookup table to find index of a given monomial in GammaMonomials (much faster than actually checking)
	% We interpret each monomial as an integer in base numOps+1 then use this integer to look in a sparse matrix for the index
	monomialIndexGammaLookup = sparse((numOps+1)^(2*maxLevelS),1);

	% Just compute for upper triangle since we can assume moment matrix is symmetric
	for i = 1:numMonomialsS
		for j = i:numMonomialsS
			GammaMonomials(:,i,j) = simplifyMonomial([flip(S(i,:)), S(j,:)], states, measurements, pure, projective, classical);
			monomialID = fromSeveralBases(GammaMonomials(:,i,j).',(numOps+1)*ones(1,2*maxLevelS)) + 1;
			% if we haven't yet seen this monomial, add to lokup table
			if monomialIndexGammaLookup(monomialID) == 0 
				monomialIndexGammaLookup(monomialID) = sub2ind([numMonomialsS, numMonomialsS],i,j);
			end
		end
	end

	% Define a lookup function that checks if a table of monomials (one per row) are contained in GammaMonomials.
	% Returns the (linear) indices.
	monomialIndexInGamma = @(monomials) monomialIndices(monomials, monomialIndexGammaLookup, numOps);

	vdisp(['Finished computing and simplifying the monomial products for moment matrix elements (', num2str(toc(t)), 's)']);

	%% Compute the index matrix

	t = tic;

	% Calculate base index matrix
	[~,~,indexMatrix] = unique(reshape(GammaMonomials, 2*maxLevelS, []).', 'stable', 'rows');
	indexMatrix = reshape(indexMatrix, numMonomialsS, numMonomialsS);

	% If we know numerical value of some terms, replace these by a common index to reduce the number of variables
	% Use numbers bigger than largest possible index to avoid clash with existing indices
	% We do this to avoid lots of constraints: if all such elements correspond to a single sdpvar, we only need 1 constraint
	ZERO = numel(indexMatrix)+1;
	ONE = numel(indexMatrix)+2;

	for i = 1:numMonomialsS
		for j = i:numMonomialsS
			monomial = GammaMonomials(:,i,j).';
			% If monomial is a state, then its trace is 1
			stateMonomials = [zeros(length(states),2*maxLevelS-1), states.'];
			if ismember(monomial, stateMonomials, 'rows')
				indexMatrix(i,j) = ONE;
			elseif projective
			% If projective hierarchy, then E_{b|y}*E{b'|y} = 0
				for k = 1:length(monomial)
					% We should check first with last element too!
					nextK = mod(k,length(monomial))+1;
					if isMeasurement(monomial(k)) && isMeasurement(monomial(nextK))
						yb1 = measurementLabels(monomial(k));
						yb2 = measurementLabels(monomial(nextK));
						if yb1(1) == yb2(1) && yb1(2) ~= yb2(2)
							indexMatrix(i,j) = ZERO;
							break
						end
					end
				end
			end
		end
	end

	% Fill in the bottom triangle (moment matrix is symmetric)
	indexMatrix = triu(indexMatrix) + triu(indexMatrix,1).';

	vdisp(['Finished computing index matrix: ', num2str(length(unique(indexMatrix))), ' unique elements (', num2str(toc(t)), 's)']);

	%% Assemble the localising matrices for the operational equivalences
	% We do this now so that, if some monomials not contained in S, we find out before trying to symmetrise

	t = tic;
	
	% First the preparation operational equivalences
	if R > 0
		% First we compute the monomial lists for the  matrices with elements 
		% tr(L(i)'*sigma(r)*L(j)) and tr(L(i)'*rho_x*L(j)) (for every state x)
		LocSigmas = zeros(2*maxLevelS, numMonomialsL, numMonomialsL, R); 
		LocRhos = zeros(2*maxLevelS, numMonomialsL, numMonomialsL, nX);
		for i = 1:numMonomialsL
			for j = i:numMonomialsL % We again only need upper triangle
				% We start with the sigma(r)
				for r = 1:R
					monomial = simplifyMonomial([flip(L(i,:)), sigmas(r), L(j,:)], states, measurements, pure, projective, classical);
					% We pad to left with identities to match length of Gamma monomials, so that we can identify terms in Gamma
					LocSigmas(:,i,j,r) = [zeros(1,2*maxLevelS-length(monomial)), monomial];
				end
				% And now the states
				for x = 1:nX
					monomial = simplifyMonomial([flip(L(i,:)), states(x), L(j,:)], states, measurements, pure, projective, classical);
					LocRhos(:,i,j,x) = [zeros(1,2*maxLevelS-length(monomial)), monomial];
				end
			end
		end

		% Identify where these monomials are in the main moment matrix, i.e. the submatrix embedding
		LocSigmaMaps = zeros(numMonomialsL, numMonomialsL, R);
		for r = 1:R
			LocSigmaMaps(:,:,r) = reshape(monomialIndexInGamma(reshape(LocSigmas(:,:,:,r), 2*maxLevelS, []).'), [numMonomialsL,numMonomialsL]);
			% Fill in the bottom triangle (localising matrix is symmetric)
			LocSigmaMaps(:,:,r) = triu(LocSigmaMaps(:,:,r)) + triu(LocSigmaMaps(:,:,r),1).';
		end
		LocRhoMaps = zeros(numMonomialsL, numMonomialsL, nX);
		for x = 1:nX
			LocRhoMaps(:,:,x) = reshape(monomialIndexInGamma(reshape(LocRhos(:,:,:,x), 2*maxLevelS, []).'), [numMonomialsL,numMonomialsL]);
			LocRhoMaps(:,:,x) = triu(LocRhoMaps(:,:,x)) + triu(LocRhoMaps(:,:,x),1).';
		end
	end

	%% Assemble the localising matrices for enforcing positivity of states and/or measurements

	if ~pure
		if isequal(L,O) && R > 0
			% Then we've already calculated these above!
			LocStatePosMaps = LocRhoMaps;
		else
			LocStatePosMaps = zeros(numMonomialsO, numMonomialsO, nX);
			for x = 1:nX
				LocState = zeros(2*maxLevelS, numMonomialsO, numMonomialsO);
				for i = 1:numMonomialsO
					for j = i:numMonomialsO
						monomial = simplifyMonomial([flip(O(i,:)), states(x), O(j,:)], states, measurements, pure, projective, classical);
						LocState(:,i,j) = [zeros(1,2*maxLevelS-length(monomial)), monomial];
					end
				end
				
				LocStatePosMaps(:,:,x) = reshape(monomialIndexInGamma(reshape(LocState, 2*maxLevelS, []).'), [numMonomialsO,numMonomialsO]);
				LocStatePosMaps(:,:,x) = triu(LocStatePosMaps(:,:,x)) + triu(LocStatePosMaps(:,:,x),1).';
			end
		end
	end

	vdisp(['Finished identifying the localising matrix elements in Gamma (', num2str(toc(t)), 's)']);

	%% Determine the dependencies between monomials due to completeness of POVMs, i.e. the constraint \sum_b E_{b|y} = \id
	% Here we record which columns of the moment matrix need to be written as combinations of other columns, and how.
	% The actual enforcing of this is done later
	
	if forceCompletenessConstraints
		t = tic;
		completenessRelations = {};
		
		% Go through each column of Gamma
		for i = 1:numMonomialsS
			monomial = S(i,:);
			% In the "standard" hierarchy, things are easier =)
			% Get locations of "last" measurement operators in monomial
			lastPOVMels = [];
			for y = 1:nY
				lastPOVMels = [lastPOVMels, find(monomial == measurements(1) + y*nB - 1)];
			end
			% If there are some, then we generate a list of combinations of other measurement operators to substitute
			if ~isempty(lastPOVMels)
				% Record that we identified a completeness relation to enforce here
				completenessRelations{end+1} = cell(1,2);
				completenessRelations{end}{1} = i;
				lastPOVMels = sort(lastPOVMels);
				n = length(lastPOVMels);
				% This monomial is to be rewritten as a linear combination of nB^n other monomials
				newMonomials = zeros(nB^n,length(monomial));
				signs = zeros(nB^n,1);
				for j = 0:nB^n - 1
					% 0 represents replace with ID, l>0 represents replace M_b by M_{b-l}
					opSubs = toSeveralBases(j, nB*ones(1,n));
					nonIDpos = lastPOVMels(opSubs > 0);
					IDpos = lastPOVMels(opSubs == 0);
					signs(j+1) = (-1)^(length(nonIDpos));
					% Then we make the substitutions
					newMonomial = monomial;
					newMonomial(nonIDpos) = newMonomial(nonIDpos) - opSubs(opSubs > 0);
					newMonomial(IDpos) = [];
					newMonomial = [zeros(1,length(IDpos)), newMonomial];
					newMonomials(j+1,:) = newMonomial;
				end
				% Find out where those monomials are and then save that correspondance
				completenessRelations{end}{2} = [monomialIndexInS(newMonomials) signs];
			end
		end
		
		disp(['Finished identifying completeness relations to be enforced (', num2str(toc(t)), 's)']);
	end
	
	%% Symmetrise the indexMatrix
	
	t = tic;

	% If we're doing some kind of symmetrisation, need to set that up
	if symmetrise > 0
		 
		stateSymGens = options.symmetryGenerators{1};
		measurementSymGens = options.symmetryGenerators{2};

		% Lookup table for each generator, specifying how each generator maps each operator
		for i = 1:length(stateSymGens)
			symGenLookup{i} = [0, stateSymGens{i}, measurementSymGens{i}+nX, sigmas];
		end

		% Compute how each monomial in S is mapped under these generators
		for i = 1:length(symGenLookup)
			SGen{i} = zeros(size(S));
			for j = 1:size(S,1)
				SGen{i}(j,:) = symGenLookup{i}(S(j,:)+1); % +1 because of we start counting operators at 0 (identity)
			end
		end

		% Determine how S has been permuted
		for i = 1:length(SGen)
			[~,symGen{i}] = ismember(S,SGen{i},'rows');
			symGen{i} = symGen{i}.';
		end

		% Symmetrise the index matrix and obtain the moment matrix sdpvar
		if symmetrise == 2 && forceCompletenessConstraints
			vdisp('No advantage from block-diagonalising when POVM completeness relations must be enforced. Reverting to standard symmetrisation');
			symmetrise = 1;
		end
		if symmetrise == 2
			% Use replab to block diagonalise
			Gamma = replab.CommutantVar.fromIndexMatrix(indexMatrix, symGen, 'symmetric', 'real');
			numSdpVars = Gamma.nbVars;
		else 
			% Or just make the indexMatrix invariant under the group and use that full symmetrised moment matrix
			uniqueIndexMatrix = replab.cvar.symmetrizeIndexMatrix(indexMatrix, symGen, 'symmetric');
			numSdpVars = max(max(uniqueIndexMatrix));
			vars = sdpvar(numSdpVars,1);
			Gamma = vars(uniqueIndexMatrix);
		end

		vdisp(['Finished creating symmetrised moment matrix: ', num2str(numSdpVars), ' variables (', num2str(toc(t)), 's)']);

	else
	%% Or construct Gamma directly without symmetrising

		% Need to make the labelling canonical
		[~,~,uniqueIndexMatrix] = unique(indexMatrix);
		uniqueIndexMatrix = reshape(uniqueIndexMatrix,size(indexMatrix));
		numSdpVars = max(max(uniqueIndexMatrix));

		vars = sdpvar(numSdpVars,1);
		Gamma = vars(uniqueIndexMatrix);

		vdisp(['Finished creating moment matrix sdpvar (', num2str(toc(t)), 's)']);

	end
	
	%% Now define the localising matrix sdpvars using the maps we identified earlier

	for r = 1:R
		for k = 1:length(E_P{r}) % for each set
			LambdaP{r}{k} = Gamma(LocSigmaMaps(:,:,r));
			for x = 1:size(E_P{r}{k},2) % and each state in that set
				LambdaP{r}{k} = LambdaP{r}{k} - E_P{r}{k}(2,x)*Gamma(LocRhoMaps(:,:,E_P{r}{k}(1,x)));
			end
		end
	end

	if ~pure
		UpsilonP = cell(1,nX);
		for x = 1:nX
			UpsilonP{x} = Gamma(LocStatePosMaps(:,:,x));
		end
	end

	%% Now set up and solve the sdp

	t = tic;

	constraints = [];

	% Constraints on main moment matrix
	constraints = [constraints, Gamma >= 0];
	% We only need to impose this for first such element, since by construction all such elements have same underlying sdpvar
	constraints = [constraints, Gamma(find(indexMatrix == ONE,1)) == 1]; 
	if projective
		% If not projective we have no such constraints!
		constraints = [constraints, Gamma(find(indexMatrix == ZERO,1)) == 0]; 
	end
	
	% Now we enforce the completeness constraints arising from \sum_b E_b = \id that we identified earlier
	% Note that if one considers projective measurements, by dropping this constraint one has a potentially weaker hierarchy but,
	% in practice, it performs just as well (the constraint is typically satisfied by the optimal Gamma anyway)
	if forceCompletenessConstraints
		tic
		for i = 1:length(completenessRelations)
			col = completenessRelations{i}{1};
			depndCols = completenessRelations{i}{2}(:,1);
			signs = completenessRelations{i}{2}(:,2);
			constraints = [constraints, Gamma(:,col) == Gamma(:,depndCols)*signs];
		end
		vdisp(['Finished applying completeness relations (', num2str(toc), 's)']);
	else 
		vdisp('Constraints arising from POVM completeness not enforced for projective measurements; set ''forceCompletenessConstraints = true'' if desired');
	end

	% Constraints on localising matrices
	for r = 1:R
		for k = 1:length(E_P{r})
			constraints = [constraints, LambdaP{r}{k} >= 0];
		end
	end
	if ~pure
		for x = 1:nX
			constraints = [constraints, UpsilonP{x} >= 0];
		end
	end
	
	% To enforce any constraints on the distribution, need to isolate p(b|x,y) terms in Gamma
	p = sdpvar(nB,nX,nY);
	for x = 1:nX
		for y = 1:nY
			for b = 1:nB
				probMonomial = simplifyMonomial([states(x), measurementFromLabels([y,b])], states, measurements, pure, projective, classical);
				pbxyIndex = monomialIndexInGamma([zeros(1,2*maxLevelS-length(probMonomial)), probMonomial]);
				p(b,x,y) = Gamma(pbxyIndex(1));
			end
		end
	end
	p = p(:);
	
	% We want to minimise tr(sigma)/K, which bounds the guessing probability
	tr_sigma = Gamma(monomialIndexInGamma([zeros(1,2*maxLevelS-1), sigmas(1)]));
	G = tr_sigma/length(E_P{1});
	
	% The constraints are given by p_constraints = {A,z}, as A*p=z
	% If p_constraints is empty, then we assume there are none.
	if ~isempty(p_constraints)
		constraints = [constraints, A*p == z];
	end

	if forceCompletenessConstraints
		vdisp(['Finished setting up SDP (',num2str(toc(t)), 's, including applying completeness relations). Starting SDP...']);
	else
		vdisp(['Finished setting up SDP (',num2str(toc(t)), 's). Starting SDP...']);
	end

	% solve the sdp
	t = tic;
	% Ensure dimacs is on!
	yalmipOptions.dimacs = 1;
	yalmipOut = optimize(constraints, G, yalmipOptions);

	if yalmipOut.problem == 0
		vdisp(['Succesfully solved SDP (', num2str(toc(t)), 's)']);
		vdisp(['Lower bound on simulation cost: ', num2str(log2(length(E_P{1})) + log2(value(G))), ' (total running time: ', num2str(toc(t_total)),'s, max dimacs measure: ', num2str(max(abs(yalmipOut.dimacs))), ')']);
	else
		vdisp(['Problems encountered solving SDP: ', yalmipOut.info, '; consult ''yalmipOut'' for more information.']);
	end
	
	G_out = value(G);
	cost = log2(length(E_P{1})) + log2(G_out);
	
	GammaOpt = value(Gamma);
	for r = 1:R
		for k = 1:length(E_P{r})
			LambdaOpt{1}{r}{k} = value(LambdaP{r}{k});
		end
	end
	if ~pure
		for x = 1:nX
			UpsilonOpt{1}{x} = value(UpsilonP{x});
		end
	else
		UpsilonOpt{1} = {};
	end
	
	function vdisp(str)
	% vdisp Displays str if the verbose flag is off
		if options.verbose
			disp(str);
		end
	end
	
end

%% Calculate and simplify the monomial list s
% Here we simplify using projectivity of measurements and/or states, so elements that reduce to same monomial
% will be identified when we generate the index matrix
function s = simplifyMonomial(s, states, measurements, pure, projective, classical)
	% Remove the zeros (identities); we'll put them back once we find the canonical form
	n = length(s);
	s = s(find(s)); 

	if classical % if classical, then everything commutes
		s = sort(s);
	end

	% If pure hierachy, then rho_x*rho_x = rho_x
	% If projective hierarchy, then M(b,y)*M(b,y) = M(b,y)
	% Warning: assumes measurements/states are contiguous list (this is ensured by construction in out implementation)
	% This is time consuming, so broken into diff loops depending on what we need to check
	i = 1;
	if pure && projective
		while length(s) > 1 && i <= length(s)
			next_i = mod(i,length(s))+1; % should check cyclicly too, since trace is cyclic
			if measurements(1) <= s(i) && s(i) <= measurements(end) && s(i) == s(next_i)
				s(i) = []; % get rid of the zeros; no need to update i
			elseif states(1) <= s(i) && s(i) <= states(end) && s(i) == s(next_i)
				s(i) = [];
			else
				i = i + 1;
			end
		end
	elseif pure
		while length(s) > 1 && i <= length(s)
			next_i = mod(i,length(s))+1;
			if states(1) <= s(i) && s(i) <= states(end) && s(i) == s(next_i)
				s(i) = [];
			else
				i = i + 1;
			end
		end	
	elseif projective
		while length(s) > 1 && i <= length(s)
			next_i = mod(i,length(s))+1;
			if measurements(1) <= s(i) && s(i) <= measurements(end) && s(i) == s(next_i)
				s(i) = [];
			else
				i = i + 1;
			end
		end		
	end

	% Convert to canonical form: amongst all cyclic permutations (and of the reverse too, since we assume a real moment matrix), 
	% choose lexicographically first
	if length(s) > 1
		sEquivClass = toeplitz(s([1 end:-1:2]), s);
		sEquivClass = [sEquivClass; toeplitz(s([end 1:end-1]), flip(s))];
		sEquivClass = sortrows(sEquivClass);
		s = sEquivClass(1,:);
	end

	% Add the zeros back to the start
	s = [zeros(1,n-length(s)), s];
end

%% Helper function
function indices = monomialIndices(monomials,monomialIndexLookup,numOps)
% monomialIndices Returns the indices of a given list of monomials (one per row) from a lookup table containing those indices
	monomialLength = size(monomials,2);
	keys = zeros(1,size(monomials,1));
	for i = 1:size(monomials,1)
		keys(i) = fromSeveralBases(monomials(i,:),(numOps+1)*ones(1,monomialLength)) + 1;
	end
	indices = full(monomialIndexLookup(keys));
	found = indices > 0;
	assert(all(found),['Error: some monomials not found in list. First such monomial: [',num2str(monomials(find(found==0,1),:)),']']);
end