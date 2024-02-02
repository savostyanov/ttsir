% Compute log-likelihoods of one rewired link in a smallworld
check_tt;
addpath('inferdata');

irun = parse_parameter('Dataset number', 1);
source_ex = parse_parameter('left end of the node after rewiring', 1);
target_ex = parse_parameter('right end of the node after rewiring', 8);
d = parse_parameter('dimension (number of nodes)', 15);

beta = 1; % 'infection rate beta'
gamma = 0.5; % 'recovery rate gamma'
delta = 0.01; % Bath infection rate
T = 10; % 'Time interval'  % !!! Specific to saved Xobs
tol = 1e-6; % 'TT approximation tolerance'

network_type = 17; % small world with 1 rewire + bath

G = WattsStrogatzDet(d,2,source_ex,target_ex);
[iadj, jadj, ~] = find(G.adjacency);
iadj = iadj';
jadj = jadj';
W_ex = sparse(iadj, jadj, ones(numel(iadj),1), d, d);
assert(norm(W_ex-W_ex',1)==0, 'W is not symmetric')
[iadj, jadj, ~] = find(W_ex);

imax_ex = adj_to_ind(W_ex);

% Initial state
n = 2*ones(d,1); % always
x0 = zeros(d,1); x0(1)=1;  % first individual infected

% Time steps of interest
tfix = (0:0.1:T)';


simulate_load_infer_data;

% Compute log-likelihoods for all possible rewirings
L = [];
for source=1:d
    for target=1:d
        G = WattsStrogatzDet(d,2,source,target);
        [iadj, jadj, ~] = find(G.adjacency);
        iadj = iadj';
        jadj = jadj';
        W = sparse(iadj, jadj, ones(numel(iadj),1), d, d);
        assert(norm(W-W',1)==0, 'W is not symmetric')
        ind = adj_to_ind(W);
        L(source,target) = cme_likelihood_si_bath(ind, beta, gamma, delta, tfix, tol, Xobs)
    end
    save(sprintf('Net%d-d%d-irun%d.mat', network_type, d, irun));
end
