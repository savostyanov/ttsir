% SIS with a bath: small probability of anyone getting infected
% Linear chain
% methods: binary mask: TT + 2 * MCMC_Shuffle + 4 * MCMC_simple
check_tt;
addpath('inferdata');

irun = parse_parameter('Dataset number', 1);
methods = parse_parameter('Method (2 = MCMC NoR, 4 = MCMC)', 2);
d = parse_parameter('dimension (number of nodes)', 9);

beta = 1; % 'infection rate beta'
gamma = 0.5; % 'recovery rate gamma'
delta = 0.01; % Bath infection rate
T = 200; % 'Time interval'  % !!! Specific to saved Xobs
tol = 1e-6; % 'TT approximation tolerance'
Nepochs = 20; % number of epochs in MCMC

network_type = 11; % chain network + bath

% Nonzero adjacency elements
iadj=[1:d-1  2:d];
jadj=[2:d    1:d-1];
XStart = 1:d;
YStart = zeros(1,d);
W_ex = spdiags(ones(d,1),1,d,d);
W_ex = W_ex + W_ex';

imax_ex = adj_to_ind(W_ex);

% Initial state
n = 2*ones(d,1); % always
x0 = zeros(d,1); x0(1)=1;  % first individual infected

% Time steps of interest
tfix = (0:0.1:T)';

simulate_load_infer_data;

infer_mcmc;
