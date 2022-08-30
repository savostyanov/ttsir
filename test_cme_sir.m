% function [outputs, ttimes_tt_full] = test_cme_sir(beta_var)
check_tt;

beta_var = [1 1 1];
d = parse_parameter('dimension', 9);
beta = parse_parameter('infection rate beta', 1);
gamma = parse_parameter('recovery rate gamma', 0.3);
T = parse_parameter('Time interval', 25);
tol = parse_parameter('TT approximation tolerance (0 to skip)', 1e-5);

network_type = parse_parameter('Network type: \n 1 - linear chain, \n 2 - 2 chains fully connected, \n 3 - 2 chains with 3 connections\n 4 - Austria\n', 1);
% Nonzero adjacency elements
switch (network_type)
    case 1
        iadj=[1:d-1  2:d];
        jadj=[2:d    1:d-1];
        XStart = 1:d;
        YStart = zeros(1,d);
%         XStart = [zeros(1,round(d/2)), ones(1,round(d/2))];
%         YStart = [d/2:-1:1, 1:d/2];
    case 2
        iadj = [1:2:d-3 3:2:d-1 2:2:d-2   4:2:d  1:2:d-1  2:2:d  ];
        jadj = [3:2:d-1 1:2:d-3 4:2:d   2:2:d-2  2:2:d    1:2:d-1];
        % 1 3 5 7   d-1
        % *-*-*-*-*-*
        % | | | | | |
        % *-*-*-*-*-*
        % 2 4 6 8   d
        XStart = reshape(repmat([1;2],1,d/2),1,[]);
        YStart = reshape(repmat(d/2:-1:1,2,1),1,[]);
    case 3
        iadj = [1:2:d-3 3:2:d-1  2:2:d-2 4:2:d  1 2 d-1 d    d/2   d/2-1];
        jadj = [3:2:d-1 1:2:d-3  4:2:d 2:2:d-2  2 1 d   d-1  d/2-1 d/2];
        % 1 3 5 7
        % *-*-*-*-*-*
        % |     |   |
        % *-*-*-*-*-*
        % 2 4 6 8   d
        XStart = reshape(repmat([1;2],1,d/2),1,[]);
        YStart = reshape(repmat(d/2:-1:1,2,1),1,[]);   
    case 4
        assert(d==9, 'Austrian graph should contain d=9 nodes')
        W = sparse(d,d);
        W(1,2) = 1; % V-T
        W(2,3) = 1; % T-Sa
        W(2,4) = 1; % T-K
        W(3,4) = 1; % Sa-K
        W(3,6) = 1; % Sa-O
        W(3,5) = 1; % Sa-St
        W(4,5) = 1; % K-St
        W(5,6) = 1; % St-O
        W(5,7) = 1; % St-N
        W(5,9) = 1; % St-B
        W(6,7) = 1; % O-N
        W(7,8) = 1; % N-W
        W(7,9) = 1; % N-B        
        % Laplacian
        W = W+W';
%         W = W - spdiags(sum(W,2),0,d,d);
%         W = W*0.5;
        [iadj, jadj, ~] = find(W);
        XStart = zeros(1,d);
        YStart = zeros(1,d);
        
        gamma = repmat(gamma, 1, d);
        beta = repmat(beta, 1,numel(iadj));
        beta( (iadj==2 & jadj==3) | (iadj==3 & jadj==2) ) = beta_var(1);
        beta( (iadj==3 & jadj==4) | (iadj==4 & jadj==3) ) = beta_var(2);
        beta( (iadj==2 & jadj==4) | (iadj==4 & jadj==2) ) = beta_var(3);
    case 5
        % Strange 4-dim network with loops
        iadj = [1 2  2 3  3 4  4 1];
        jadj = [2 1  3 2  4 3  1 4];
        XStart = zeros(1,d);
        YStart = zeros(1,d); 
    case 6
        % quasi-2d quasi-full layers
        ijadj = [];
        rng(0);
        p = 0.5; % probability of sampling an edge
        for i=1:d-1
            for j=1:4
                for k=1:4
                    if (rand<p)
                        ijadj = [ijadj, [j+(i-1)*4, k+(i+1-1)*4; k+(i+1-1)*4, j+(i-1)*4]];
                    end
                end
            end
        end
        ijadj = unique(ijadj', 'rows', 'stable');
        iadj = ijadj(:,1)';
        jadj = ijadj(:,2)';
        XStart = reshape(repmat((1:d), 4, 1), 1, []);
        YStart = reshape(repmat((1:4)', 1, d), 1, []);
        
        d = d*4;
end

%% Visualise the network
G = sparse(iadj, jadj, beta, d, d);
G = graph(G);
figure(2); 
h = plot(G, 'Layout', 'force', 'XStart', XStart, 'YStart', YStart, 'Iterations', 100); % , 'EdgeLabel',G.Edges.Weight);
% h.LineWidth = 3;
% h.MarkerSize = 8;
% h.NodeFontSize = 20;

h.Parent.Box=false;
h.Parent.XAxis.Visible=false;
h.Parent.YAxis.Visible=false;

% print(sprintf('graph-p%d.pdf', round(swprob*100)), '-dpdf')

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
Icritical = parse_parameter('Quantities of interest: \n 0 - mean and variance, \n Icritical>0 - probability of I>=Icritical\n', 0);
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
tfix = (0:0.1:T)';

%% TT
if (tol>0)
    % Run this to solve in TT
    ttimes_tt_full = tic;
    [Uglobal,tglobal,outputs_raw,norm_u,ttimes_tt] = tt_cme_sir(x0, iadj, jadj, beta, gamma, T, tol, 'outfun', Icritical, 'tfix', tfix, 'verbose', 0, 'max_nt', 8);
%     [Uglobal,tglobal,outputs2,ttimes_tt,XY] = tt_cme(z,wfun,x0,T, n, tol, 'outfun', outfun, 'tfix', tfix, 'verbose', 0, 'kickrank_w', 16);
    ttimes_tt_full = toc(ttimes_tt_full);
    
    ttranks = cellfun(@(x)size(x,4), Uglobal);
    
    outputs = outputs_raw ./ norm_u;
    
%     err_tt = sqrt(sum((outputs - outputs2).^2,1)./sum(outputs2.^2,1))
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
xlabel('t');
title('output statistics')

% end
