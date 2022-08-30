global swprob

swprob = parse_parameter('probability of rewiring', 0.03);
d = parse_parameter('dimension', 50);
beta = parse_parameter('infection rate beta', 1);
gamma = parse_parameter('recovery rate gamma', 0.3);
T = parse_parameter('Time interval', 25);
tol = parse_parameter('TT approximation tolerance (0 to skip)', 1e-7);

% Small world
rng(46); % interesting values 0,9,48,   6,15,    4,14,17,      17, 33, 44, 46
         % dimension           100       80        60          50
G = WattsStrogatz(d,2,swprob);
[iadj, jadj, ~] = find(G.adjacency);
iadj = iadj';
jadj = jadj';
XStart = cos((1:d)*2*pi/d);
YStart = sin((1:d)*2*pi/d);

%% Visualise the network
G = sparse(iadj, jadj, beta, d, d);
G = graph(G);
figure(2); 
h = plot(G, 'Layout', 'force', 'XStart', XStart, 'YStart', YStart, 'Iterations', 0);

h.Parent.Box=false;
h.Parent.XAxis.Visible=false;
h.Parent.YAxis.Visible=false;

print(sprintf('graph-p%d.pdf', round(swprob*100)), '-dpdf')

%% Initial state
n = 3*ones(d,1); % always
x0 = zeros(d,1); x0(1)=1;  % first individual infected

%% Populate z and w from adjacency matrix
% z is d x M
z = zeros(d, d+numel(iadj));
wfun = cell(2, d+numel(iadj));
for m=1:numel(iadj)
    % infection iadj -> jadj
    z(jadj(m),m) = 1; % increment the state by 1: only 0->1 transition is possible, exactly what we need
    wfun{1,m} = @(x)beta(min(m,numel(beta))) * double(x(:,1)==0).*double(x(:,2)==1);
    wfun{2,m} = [jadj(m), iadj(m)];
end
% recovery
for m=1:d
    z(m,numel(iadj)+m) = 1; % increment the state by 1: only 1->2 transition is possible, exactly what we need
    wfun{1,numel(iadj)+m} = @(x)gamma(min(m,numel(gamma))) * double(x==1);
    wfun{2,numel(iadj)+m} = m;
end


% Functions of quantities of interest
Icritical = parse_parameter('Quantities of interest: \n 0 - mean and variance, \n Icritical>0 - probability of I>=Icritical\n', d*4/5);
if (Icritical==0)
    % Mean and second moment of infected and recovered
    outfun = {
        @(x)sum(double(x==1),2), @(x)sum(double(x==1),2).^2, ...
%         @(x)sum(double(x==2),2), @(x)sum(double(x==2),2).^2 ...
        };    
elseif (Icritical==-1)
    % Mean number of susceptible
    outfun = {@(x)sum(double(x==0),2), @(x)sum(double(x==0),2).^2};
elseif (Icritical==-2)
    % Probability of x_{d/2} = s
    outfun = {@(x)double(x(:,d/2)==0)};
else
    % Probability of I>=Icritical
    outfun = {
        @(x)double(sum(x==1,2)>=Icritical), @(x)double(sum(x==1,2)>=(Icritical+2))
        };
end

% Time steps of interest
tfix = (0:0.01:T)';

%% TT
if (tol>0)
    % Run this to solve in TT
    ttimes_tt_full = tic;
    [Uglobal,tglobal,outputs_raw,norm_u,ttimes_tt] = tt_cme_sir(x0, iadj, jadj, beta, gamma, T, tol, 'outfun', Icritical, 'tfix', tfix, 'verbose', 0, 'max_nt', 0, 'time_error_damp', 0, 'nswp', 20, 'resid_damp', 1000, 'local_iters', 1000);
    ttimes_tt_full = toc(ttimes_tt_full);
    
    ttranks = cellfun(@(x)size(x,4), Uglobal);
    
    outputs = outputs_raw ./ norm_u;
end


%% SSA
Nsim = parse_parameter('Number of samples for SSA (0 to skip)', (~(tol>0))*1e4);
if (Nsim>0)
    % Run this for SSA
    ttimes_ssa = tic;
    outputs = ssa(z,wfun,x0,T,Nsim,tfix, 100, outfun);
    ttimes_ssa = toc(ttimes_ssa);
%     outputs = expectations_ssa(outputs, outfun);
end

%% Postprocess and plot stuff
if (Icritical<=0)
    % compute std dev
    outputs(:,2:2:end) = sqrt(abs(outputs(:,2:2:end) - outputs(:,1:2:end-1).^2));
end

figure(1); plot(tfix, outputs);

