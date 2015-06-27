%%Simulate single graph
nrepli = 150; %number of replicates
p = 100; %number of nodes
prob = 0.1; %sparsity level
opts = []; 
opts.maxv = 0.4; %maximum value (absolute value) of the entries in the precision matrix
opts.minv = 0.1; %minimum value (absolute value) of the entries in the precision matrix
opts.add = 0.5; %the added value to the diagonal in the precision matrix
sim = simsig(nrepli, p, prob, opts);
sim = simsig(nrepli, p, prob); %opts can be omitted to use the default setting
sim.data %n by p data matrix
sim.omega %the precision matrix
sim.sup %support of the precision matrix
%%
%%Run the Bayesian Neigborhood Selection (BNS) method on single graph
opts = [];
opts.niter = 20000; %total number of iterations
opts.br = 10000; %burn-in
opts.nu = 0; %the inverse gamma prior on sigma
opts.lambda = 0; %the inverse gamma prior on sigma
opts.parallel = 0; %Use parallel computing to update each line of beta? 0: no; 1: yes.
%see manuscript for how to specify q, l and delta
q = 0.1;
l = 0.1;
delta = 0.1;
result1 = getBNS( sim.data, q, l, delta, opts );
result1.tau1 %tau_1
result1.tau0 %tau_0
result1.postprob %posterior probability of the edges

%parallel computing
opts = [];
opts.niter = 20000; %total number of iterations
opts.br = 10000; %burn-in
opts.nu = 0; %the inverse gamma prior on sigma
opts.lambda = 0; %the inverse gamma prior on sigma
opts.parallel = 1;
%Use parpool(# of cores) to enable parallele computing
parpool(4); %4 cores are used
%see manuscript for how to specify q, l and delta
q = 0.1;
l = 0.1;
delta = 0.1;
result2 = getBNS( sim.data, q, l, delta, opts );
%%
%%Simulate multiple graphs
nrepli = 150; %number of replicates
p = 100; %number of nodes
ngraph = 3; %number of graphs
prob = 0.1; %sparsity level
change = 0.2; %1-the level of sharing, 0 means the graphs share all edges, 1 means the graphs share no more edges than random
value = 'different'; %'different' or 'same': whether shared edges take the same value in the precision matrix
opts = []; 
opts.maxv = 0.4; %maximum value (absolute value) of the entries in the precision matrix
opts.minv = 0.1; %minimum value (absolute value) of the entries in the precision matrix
opts.add = 0.5; %the added value to the diagonal in the precision matrix
sim = simS(nrepli, p, ngraph, prob, change, value, opts);
sim = simS(nrepli, p, ngraph, prob, change, value); %opts can be omitted to use the default setting
sim.dataA %data matrix
sim.supA %support of the precision matrix
sim.omega %precision matrix 

%%
%%Run the Bayesian Neigborhood Selection (BNS) method on multiple graphs
opts = [];
opts.niter = 20000; %total number of iterations
opts.br = 10000; %burn-in
opts.nu = 0; %the inverse gamma prior on sigma
opts.lambda = 0; %the inverse gamma prior on sigma
opts.eta1 = -0.5; %the value for eta_1
opts.etaS = 0; %the initial value of eta_S
opts.fixetaS = 0; %fix or update eta_s by Metropolis-Hastings algorithm
opts.parallel_graph = 0; %Use parallel computing when updating multiple graphs? 0: no; 1: yes.
opts.parallel_line = 0; %Use parallel computing when updating each line of beta? 0: no; 1: yes.
%important: opts.parallel_graph and opts.parallel_line cannot both be 1
%see manuscript for how to specify l and delta
l = 0.1;
delta = 0.1;
gA = -1; %the initial latent states to sample from, -1 means random
result = getBNSspatial( sim.dataA, l, delta, gA, opts);
result.tau1 %tau_1
result.tau0 %tau_0
result.postprob %posterior probability of the edges
result.etaSA %trace of etaS
%%
%%Simulate multiple graphs with temporal dependency, the graph structure
%%evolves over time by hidden markov model
nrepli = 100; %number of replicates
p = 50; %number of nodes
nt = 10; %number of time points
prob = 0.1; %level of sparsity
tran = 0.2; %transition probability
value = 'different'; %'different' or 'same': whether shared edges take the same value in the precision matrix
opts = []; 
opts.maxv = 0.4; %maximum value (absolute value) of the entries in the precision matrix
opts.minv = 0.1; %minimum value (absolute value) of the entries in the precision matrix
opts.add = 0.5; %the added value to the diagonal in the precision matrix
sim = simT(nrepli, p, nt, prob, tran, value, opts);
sim.dataA %data matrix
sim.supA %support of the precision matrix
sim.omega %precision matrix 

%%
%%Run the Bayesian Neigborhood Selection (BNS) method on multiple graphs
%%with temporal dependency
opts = [];
opts.niter = 20000; %total number of iterations
opts.br = 10000; %burn-in
opts.nu = 0; %the inverse gamma prior on sigma
opts.lambda = 0; %the inverse gamma prior on sigma
opts.eta1 = -0.5; %the value for eta_1
opts.etaT = 0; %the initial value of eta_T
opts.fixetaT = 0; %fix or update eta_T by Metropolis-Hastings algorithm
opts.parallel_graph = 0; %Use parallel computing when updating multiple graphs? 0: no; 1: yes.
opts.parallel_line = 0; %Use parallel computing when updating each line of beta? 0: no; 1: yes.
%important: opts.parallel_graph and opts.parallel_line cannot both be 1
%see manuscript for how to specify l and delta
l = 0.1;
delta = 0.1;
gA = -1; %the initial latent states to sample from, -1 means random
result = getBNStemporal( sim.dataA, l, delta, gA, opts);
result.tau1 %tau_1
result.tau0 %tau_0
result.postprob %posterior probability of the edges
result.etaSA %trace of etaS

%%
%%Simulate multiple graphs with temporal and spatial dependency 
nrepli = 100; %number of replicates
p = 50; %number of nodes
nt = 10; %number of time points
nbr = 3; %number of loci
prob = 0.1; %level of sparsity
tran = 0.4; %transition probability
purt = 0.2; %purturbation
value = 'different'; %'different' or 'same': whether shared edges take the same value in the precision matrix
opts = []; 
opts.maxv = 0.4; %maximum value (absolute value) of the entries in the precision matrix
opts.minv = 0.1; %minimum value (absolute value) of the entries in the precision matrix
opts.add = 0.5; %the added value to the diagonal in the precision matrix
sim = simST(nrepli, p, nt, nbr, prob, tran, purt, value, opts);
sim.dataA %data matrix
sim.supA %support of the precision matrix
sim.omega %precision matrix 

%%
%%Run the Bayesian Neigborhood Selection (BNS) method on multiple graphs
%%with temporal and spatial dependency 
opts = [];
opts.niter = 20000; %total number of iterations
opts.br = 10000; %burn-in
opts.nu = 0; %the inverse gamma prior on sigma
opts.lambda = 0; %the inverse gamma prior on sigma
opts.eta1 = -0.5; %the value for eta_1
opts.etaS = 0; %the initial value of eta_S
opts.etaT = 0; %the initial value of eta_T
opts.fixetaT = 0; %fix or update eta_T by Metropolis-Hastings algorithm
opts.fixetaS = 0; %fix or update eta_S by Metropolis-Hastings algorithm
opts.parallel_t = 0; %Use parallel computing when updating graphs over time? 0: no; 1: yes.
opts.parallel_s = 0; %Use parallel computing when updating graphs over space? 0: no; 1: yes.
opts.parallel_line = 0; %Use parallel computing when updating each line of beta? 0: no; 1: yes.
%important: only one of opts.parallel_t, opts.parallel_s and opts.parallel_line can be 1
%see manuscript for how to specify l and delta
l = 0.1;
delta = 0.1;
gA = -1; %the initial latent states to sample from, -1 means random
result = getBNSst( sim.dataA, l, delta, gA, opts);
result.tau1 %tau_1
result.tau0 %tau_0
result.postprob %posterior probability of the edges
result.etaSA %trace of etaS
result.etaTA %trace of etaT
