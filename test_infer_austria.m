% SIS with a bath: small probability of anyone getting infected
% Austria network
% methods: binary mask: TT + 2 * MCMC_Shuffle + 4 * MCMC_simple
check_tt;
addpath('inferdata');

irun = parse_parameter('Dataset number', 1);
methods = parse_parameter('Method (2 = MCMC NoR, 4 = MCMC)', 2);
d = parse_parameter('dimension (number of nodes)', 9);

beta = 1; % 'infection rate beta'
gamma = 0.5; % 'recovery rate gamma'
delta = 0.01; % Bath infection rate
T = 1000; % 'Time interval'  % !!! Specific to saved Xobs
tol = 1e-6; % 'TT approximation tolerance'
Nepochs = 20; % number of epochs in MCMC

network_type = 14; % austria + bath

% Nonzero adjacency elements
assert(d==9, 'Austrian graph should contain d=9 nodes')
W_ex = sparse(d,d);
W_ex(1,2) = 1; % V-T
W_ex(2,3) = 1; % T-Sa
W_ex(2,4) = 1; % T-K
W_ex(3,4) = 1; % Sa-K
W_ex(3,6) = 1; % Sa-O
W_ex(3,5) = 1; % Sa-St
W_ex(4,5) = 1; % K-St
W_ex(5,6) = 1; % St-O
W_ex(5,7) = 1; % St-N
W_ex(5,9) = 1; % St-B
W_ex(6,7) = 1; % O-N
W_ex(7,8) = 1; % N-W
W_ex(7,9) = 1; % N-B
% Laplacian
W_ex = W_ex+W_ex';
[iadj, jadj, ~] = find(W_ex);

imax_ex = adj_to_ind(W_ex);

% Initial state
n = 2*ones(d,1); % always
x0 = zeros(d,1); x0(1)=1;  % first individual infected

% Time steps of interest
tfix = (0:0.1:T)';

simulate_load_infer_data;

infer_mcmc;
