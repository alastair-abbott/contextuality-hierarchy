# Bounding and simulating contextual correlations in quantum theory

**Authors: Armin Tavakoli, Emmanuel Zambrini Cruzeiro, Roope Uola, and Alastair A. Abbott**

This code provides an implementation of a SDP hierarchy of necessary conditions to bound the set of quantum correlations in contextuality scenarios, as described in

> [A. Tavakoli, E. Zambrini Cruzeiro, R. Uola, A. A. Abbott, "Bounding and simulating contextual correlations in quantum theory", PRX Quantum **2**, 020334 (2021)](https://doi.org/10.1103/PRXQuantum.2.020334)

and can handle almost completely arbitrary contextuality scenarios, with many options, for example to restrict to quantum models comprised of pure states.
A similar hierarchy to bound the simulation cost of quantum or post-quantum contextual correlations is also provided.

## Installation

This code requires the following:
- [MATLAB](https://fr.mathworks.com/products/matlab.html) (tested with R2018a and newer; probably compatible with older releases too, but not tested);
- [Yalmip](https://yalmip.github.io/) (available in your matlab path);
- An sdp solver such as [MOSEK](https://www.mosek.com/), [SeDuMi](https://github.com/sqlp/sedumi), or [SCS](https://github.com/bodono/scs-matlab), installed and working with Yalmip.

In addition, if you want to exploit symmetries in your problem (which greatly reduces the size of the SDP problem), one must install:
- [RepLAB](https://replab.github.io/replab/) (with the replab folder in your matlab path) 
    - Note: a recent developer build is required; the release 0.9.0 will not work. We recommend cloning the replab git repository following the advanced instructions [here](https://replab.github.io/replab/tutorials/installation.html) (Option 2). We have tested commit [6f2fdd](https://github.com/replab/replab/tree/6f2fdda68e4e36ee2adfa28fd88c520f166b5040).

The contextuality hierarchy can simply be run from wherever you wish to install it. The examples can be run as is (as long as the base directory is in added to your path), and provide the easiest way to examine the usage of the SDP hierarchy.

## Usage

The main functionality of the contextuality hierarchy is provided by the function [`boundContextualCorrelations`](boundContextualCorrelations.m), and its usage is probably best discerned by looking at the included [examples](examples). Here a brief summary of the required input and options is provided. The notation and terminology follows the paper cited at the start of this readme.

The function header for running the hierarchy is `[objOpt, yalmipOut, GammaOpt, LambdaOpt, UpsilonOpt] = boundContextualCorrelations(n_xyb, E_P, E_M, p_constraints, w, levels, options)`.

Firstly, one must provide a specification of the "contextuality scenario" being considered. This consists of:
- `n_xyb = [nX, nY, nB]`, specifying the number of preparations, measurements, and measurement outcomes, respectively 
- `E_P,E_M`, the preparation and measurement operational equivalences
    - Following the paper, `E_P` is a cell array of operational equivalences `E_P{r}`, which itself is a $`1\times K_r`$ dimensional cell array so that `E_P{r}{k}` is a 2-row matrix. The first specifies the preparations (numbered `1` to `nX`) in the set $`S^{(r)}_k`$, and the second sepcifies the corresponding weights $`\xi^{(r)}_k`$.
    - `E_M` is defined specified analagously
    - Both these arguments must be specified, but can be empty (e.g. `{}`) if no operational equivalences of the corresponding type are present in the contextuality scenario.

Constraints that the observed probability distribution $`p(b|x,y)`$ should satisfy can be specified in the argument `p_constraints`:
- If no particular constraints are to be enforced (as is the case if one simply wishes to maximise the score for a non-contextuality inequality, for example), then one can simply set `p_constraints = {}`.
- Otherwise, one should pass `p_constraints = {A,z}`, where `A` is a $`m\times (n_X n_Y n_B)`$ matrix and `z` is an $`m`$ element column vector such that the probability vector $`p=(p(1|1,1),\dots,p(n_B|1,1),p(1|2,1),\dots)^T`$ satisfies $`A p = z`$.
    - This allows one, for example, to fix completely the probability distribution and ask the feasibility problem of whether such a distribution is compatible with a quantum model.

The objective function can be specified either as a linear functional of the probability distribution, or instead a robust feasibility problem can be solved:
- One can maximise the function $`\sum_{b,x,y}w(b,x,y)p(b|x,y)`$ by passing the $`n_B\times n_X \times n_Y`$ element array `w` with elements `w(b,x,y)`.
- To solve a feasibility problem, pass `w = []` or `w = zeros(nB,nX,nY)`.
    - In this case, the hierarchy will solve the more robust feasibility problem by maximising $`\lambda`$ such that $`\Gamma - \lambda \bm{1} \succeq 0`$, $`\Upsilon^x - \lambda \bm{1} \succeq 0`$, and $`\Upsilon^{(b,y)} - \lambda \bm{1} \succeq 0`$. A strictly negative solution (i.e., beyond numerical precision) hence provides a measure of the infeasibility of the problem.

One must specify the levels of the SDP hierarchy; i.e., what moment lists $`\mathcal{S},\mathcal{L},\mathcal{O}`$ one wishes to use. These are specified in a high-level manner by signalling what "combinations" of different operators are to be combined. 
- We use the notation that `1` stands in for all states, `2` for all measurements, and `3` for the operators $`\sigma_r,\tau_q`$ used to enforce the operational equivalences.
- A level specification is thus an array over $`\{1,2,3\}`$. For example, `[1 2 1 2]` means "all monomials of the form $`\rho_x E_{b|y}\rho_{x'}E_{b'|y'}`$", for any $`x,x',y,y',b,b'`$.
- One thus specifies `levels = {levelsS, levelsL, levelsO}`, where each of of these is a cell array of level specifications. E.g., one may have `levelsS = {1,2,3,[1 1],[1 3],[2 2], [2 3]}`, which defines a sublevel of the "2nd level" of the hierarchy.
    - Note: one must make sure that the elements of the localising matrices are elements of the moment matrix. This notably means that `levelsS` must contain larger level specifications that `levelsL` and `levelsO`. The code will detect if this is not the case and throw an error, indicating the monomials that were not found.
    - Internally, the code labels the operators so that the identity is `0`, the states are `1` to `nX`, the measurements are `nX+1` to `nX+nY*nB`, and then one has the $`\sigma_r`$ and $`\tau_q`$. This may be helpful to identify what monomials could not be found.

The hierarchy also takes several optional options, which can be provided in a structure as the final, optional, argument.
- `verbose` (default: `false`): Print out helpful information and progress as the code runs.
- `classical` (default: `false`): Treat all operators as commuting and thus classical.
- `pure` (default: `false`): Treat all states as being pure.
- `projective` (default: `false`): Treat all measurements as projective.
- `yalmipOptions`: A yalmip `sdpsettings()` structure allowing the solver to be changed and settings to be specified (e.g., by setting the solver to verbose mode).
- `forceCompletenessConstraints` (default: `false`): If the projective version of the hierarchy is run (as will always be the case if `E_M` is empty), we don't, by default, impose the completeness of the POVMs, as in general this seems not to help. This setting allows them to nonetheless be imposed if desired.
- `sqrtRhoHierarchy` (default: `false`): Use the variant of the hierachy described in Appendix B of the paper, in which the state operators are treated as $`\sqrt{\rho_x}`$ instead of $`\rho_x`$ (e.g., in the specification of the levels). The constraints are then adjusted internally as required.
- `sqrtEHierarchy` (default: `false`): The analagous variant in which the measurement operators are interpreted as $`\sqrt{E_{b|y}}`$.
- `symmetrise` (default: `0`): Specify what symmetrisation to use.
    - `0`: No symmetrisation applied.
    - `1`: Make the moment matrix invariant under the supplied symmetries (see below). Note that RepLAB is required for this.
    - `2`: Block diagonalise the moment matrix. This option invokes longer pre-processing but results in a significantly smaller SDP if the completeness relations are not applied. This option should be treated as experimental. Requires RepLAB.
- `symmetryGenerators = {stateSymGens, measurementSymGens}`: This option is required if `symmetrise > 0`. `stateSymGens` and `measurementSymGens` are cell arrays of the same length specifying the permutations on the states and measurements under which the objective function and constraints are invariant. These permutations should be understood as being applied simultaneously: the $`i`$th state permutation is applied at the same time as the $`i`$th measurement permutation. Each element of `stateSymGens` should thus be a permutation of $`1,\dots,n_X`$, and each element of `measurementSymGens` a permutation of $`1,\dots,n_Y*n_B`$. The measurements are understood here has being listed in the order $`E_{1|1},E_{2|1},\dots,E_{n_B|1},E_{1|2},\dots,E_{n_B|n_Y}`$.

The hierarchy function returns the array `[objOpt, yalmipOut, GammaOpt, LambdaOpt, UpsilonOpt]`.
- `objOpt`: The bound on the objective function or, if a feasibility problem is run, a measure of the feasibility of the problem.
- `yalmipOut`: Contains information useful for debugging, and notably the `dimacs` error metrics.
- `GammaOpt, LambdaOpt, UpsilonOpt`: The moment and localising matrices corresponding to the solution of the SDP.


Finally, a similar hierarchy bounding the simulation cost of contextual correlations (see paper) is provided as the function `boundSimulationCost(n_xyb, E_P, p_constraints, levels, options)`. The usage here is largely the same, except that only contextuality scenarios with one preparation operational equivalence, and no measurement operational equivalences, are supported. Similarly, no objective function is provided, as the hierarchy instead provides a lower bound on the simulation cost to obtain correlations satisfying the operational equivalence and the constraints given by `p_constraints`. The first argument returned by the function is the simulation cost in bits.


