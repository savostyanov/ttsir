function [Uglobal,tglobal,outputs,norm_u,ttimes_solve]=tt_cme_sir(x0, iadj, jadj, beta, gamma, T, tol, varargin)
% Construct and solve the SIR CME using tamen.
% Inputs: 
%   x0: deterministic initial state (d x 1 vector)
%   iadj: a vector of row indices of nonzero adjacency elements
%   jadj: a vector of column indices of nonzero adjacency elements
%   beta: infection rate
%   gamma: recovery rate
%   T: final time
%   tol: TT truncation tolerance
% Optional variables can be given in pairs 'ParameterName', ParameterValue
% Most crucial are:
%   outfun: specifies the outputs. Set outfun to 0 to compute the mean and
%       second moment of the total number of infected individuals, a positive
%       value means the threshold I_* for the probability of I>=I_*, and a
%       value of -1 will compute the mean number of susceptible individuals
%   tfix: a vector of time points where the outputs are sought.
% Other parameters are passed to tamen and amen_solve. A special treatment
% is applied to max_nt: a positive value is used as the number of Chebyshev
% grid points on an individual time interval for tamen, which is run
% adapting the time intervals; whereas max_nt=0 invokes an implicit Euler
% discretization on the time steps given in tfix instead.
% Outputs:
%   Uglobal: a cell array of TT formats of the solution from tamen
%   tglobal: a cell array of time steps from tamen
%   outputs: a matrix of statistics (n_time_steps x n_outputs)
%   norm_u: a vector of normalising constants (n_time_steps x 1)
%   ttimes_solve: CPU time of the solution
global swprob

atol = 1e-14;
opts.max_nt = 8;
opts.kickrank = 4; % for tamen
opts.nswp = 10; % for tamen -- should be enough

for i=1:2:length(varargin)
    switch lower(varargin{i})
        case 'atol'
            atol = varargin{i+1};
            fprintf('Found parameter: %s=%3.3e\n', varargin{i}, varargin{i+1});
        case 'outfun'
            outfun = varargin{i+1};       
            fprintf('Found parameter: %s\n', varargin{i});
        case 'tfix'
            tfix = varargin{i+1};       
            fprintf('Found parameter: %s\n', varargin{i});
            if (max(tfix)>T)
                warning('tfix>T detected, the output may be inaccurate');
            end
        otherwise
            warning('Option %s was not recognized by tt_cme. Passing to amens', varargin{i});
            opts.(varargin{i}) = varargin{i+1};
    end
end


d = numel(x0);

% Det Initial state
u = tt_unit(3,1,x0(1)+1);
for i=2:d
    u = tkron(u, tt_unit(3,1,x0(i)+1));
end

if (numel(gamma)==1)
    gamma = repmat(gamma, 1, d);
end
if (numel(beta)==1)
    beta = repmat(beta, 1, numel(iadj));
end

% Matrix parts
A = cell(1,0);
% Recovery
for i=1:d
    Ai = ([0,0,0;1,0,0;0,1,0] - eye(3)) * diag([0;1;0]) * gamma(i);
    Ai = tt_matrix(Ai);
    if (i>1)
        Ai = tkron(tt_eye(3,i-1), Ai);
    end
    if (i<d)
        Ai = tkron(Ai, tt_eye(3,d-i));
    end
%     A(:,end+1) = core2cell(Ai);
    A{1,end+1} = tt_tensor(Ai);
end
% Infection
for i=1:numel(iadj)
    % iadj -> jadj
    Ai = beta(i)*diag([0;1;0]);
    Aj = ([0,0,0;1,0,0;0,1,0] - eye(3)) * diag([1;0;0]);
    Ai = tt_matrix(Ai);
    Aj = tt_matrix(Aj);
    if (iadj(i)<jadj(i))
        if (iadj(i)>1)
            Ai = tkron(tt_eye(3,iadj(i)-1), Ai);
        end
        if (jadj(i)<d)
            Aj = tkron(Aj, tt_eye(3, d-jadj(i)));
        end
        if (jadj(i)-iadj(i)>1)
            Ai = tkron(Ai, tt_eye(3, jadj(i)-iadj(i)-1));
        end
        Ai = tkron(Ai, Aj);
    else
        if (jadj(i)>1)
            Aj = tkron(tt_eye(3,jadj(i)-1), Aj);
        end
        if (iadj(i)<d)
            Ai = tkron(Ai, tt_eye(3, d-iadj(i)));
        end        
        if (iadj(i)-jadj(i)>1)
            Aj = tkron(Aj, tt_eye(3, iadj(i)-jadj(i)-1));
        end        
        Ai = tkron(Aj, Ai);
    end
%     A(:,end+1) = core2cell(Ai);
    A{1,end+1} = tt_tensor(Ai);
end

A = amen_sum(A, ones(numel(A),1), atol, 'y0', tt_tensor(tt_eye(3,d)), 'nswp', opts.nswp_w, 'kickrank', opts.kickrank_w);
A = tt_matrix(A);
A = core2cell(A);


% Observation of total number of infected
Itot = cell(d,1);
Itot{1} = reshape([[0;1;0], [1;1;1]], 1, 3, 2);
for i=2:d-1
    Itot{i}(1,:,1) = [1 1 1];
    Itot{i}(2,:,2) = [1 1 1];
    Itot{i}(2,:,1) = [0 1 0];
    Itot{i}(1,:,2) = [0 0 0];
end
Itot{d} = [[1 1 1]; [0 1 0]];
Itot = cell2core(tt_tensor, Itot);

P0 = ItotIndicator(d, outfun);
P2 = ItotIndicator(d, outfun+2);

% Total number of susceptible
Stot = cell(d,1);
Stot{1} = reshape([[1;0;0], [1;1;1]], 1, 3, 2);
for i=2:d-1
    Stot{i}(1,:,1) = [1 1 1];
    Stot{i}(2,:,2) = [1 1 1];
    Stot{i}(2,:,1) = [1 0 0];
    Stot{i}(1,:,2) = [0 0 0];
end
Stot{d} = [[1 1 1]; [1 0 0]];
Stot = cell2core(tt_tensor, Stot);


ons = tt_ones(3,d);

outs = [];
if (~isempty(outfun))
    if isscalar(outfun) && (outfun==0)
        outs = {Itot, Itot.^2};
    end
    if isscalar(outfun) && (outfun>0)
        % Indicator of I>=outfun
        outs = {P0, P2};
    end
    if isscalar(outfun) && (outfun==-1)
        outs = {Stot};  % just number of susceptible
    end
end
Nout = numel(outs);




% Space-time solution via adaptive tamen
Uglobal = tkron(u, tt_ones(opts.max_nt));
Uglobal = tt_matrix(Uglobal, Uglobal.n, 1);
Uglobal = core2cell(Uglobal);

tref = [];
if (~isempty(tfix))
    tref = sort(tfix, 'ascend');
    opts.correct_2_norm = false;       
end

B = [A, core2cell(tt_eye(3,d))];

z=[];

tic_tt_cme = tic;
if (opts.max_nt>0)
    % Scale the matrix to desired time
    A(1,:) = cellfun(@(x)(x*T), A(1,:), 'UniformOutput', false);    
    [Uglobal,tglobal] = tamen(Uglobal,A,tol,opts);
    % Scale tglobal back to T
    tglobal = cellfun(@(t)t*T, tglobal, 'UniformOutput', false);
else
    opts.obs = []; core2cell(tt_matrix(P0,3,1));
    opts.kickrank = 8;
    % opts.obs = {Stot};
        
    ttranks = zeros(numel(tref), 1);
    norm_u = zeros(numel(tref), 1);
    norm_u(1) = dot(ons,u);
    for m=1:numel(tref)-1
        B{1,1} = -(tref(m+1)-tref(m)) * A{1,1};
%         B = A;
%         B(1,:) = cellfun(@(x)(-x*(tref(m+1)-tref(m))), A(1,:), 'UniformOutput', false);    
%         B = [B, Icell];
%         B = cell2core(tt_matrix, B);
        [u,opts,swp,  z] = amen_solve(B, u, tol, opts, u, opts.obs, z);
        tglobal{m,1} = 1;
        
%         [Uglobal(:,m),tglobal(m,1),~,u] = tamen(Uglobal(:,end),B,tol,opts);

%         [Uglobal,tglobal(m,1),~,u, z] = tamen(Uglobal,B,tol,opts, z);
%         u = tt_tensor(cell2core(tt_matrix, u));
        
        tglobal{m,1} = tref(m) + tglobal{m,1} * (tref(m+1)-tref(m));
        
        ttranks(m+1) = max(u.r);
        
        norm_u(m+1) = dot(ons,u);    
        for j=1:Nout
            outputs(m+1,j) = dot(outs{j}, u); % /norm_u(m+1); 
        end
        figure(1); plot(tfix(2:m+1), outputs(2:m+1,:));

        opts.nswp = 1;        
        
        ttimes_solve = toc(tic_tt_cme);
%         if (abs(tref(m+1)*10 - round(tref(m+1)*10))<1e-8)
%             % Save every 0.1 time units
%             save(sprintf('../scratch/ttsir/sw-N%d-p%d-t%g.mat', d, round(swprob*100), tref(m+1)), 'u', 'ttimes_solve');
%         end

        fprintf('swprob=%d, time step %d[%g] done, 1-|u|=%3.3e, rank=%d, P0=%g\n', round(swprob*100), m, tref(m+1), 1-norm_u(m+1), ttranks(m+1), outputs(m+1,1));
    end
end
ttimes_solve = toc(tic_tt_cme);


if (opts.max_nt>0)
    % Determine tfix
    if (isempty(tfix))
        tfix = cell2mat(tglobal);
    end
    outputs = zeros(numel(tfix), Nout);
    norm_u = zeros(numel(tfix), 1);
    % Iterate over all time step
    for i=1:numel(tfix)
        u = extract_snapshot(Uglobal, tglobal, tfix(i));
        u = cell2core(tt_matrix, u);
        u = tt_tensor(u);
        norm_u(i) = dot(ons,u);
        for j=1:Nout
            outputs(i,j) = dot(outs{j}, u); % /norm_u;
        end
    end
    
    fprintf('t=%g, 1-|u|=%3.3e, outputs:\n\t', tfix(i), 1-norm_u(i));
    for j=1:Nout
        fprintf('%3.5e  ', outputs(i,j));
    end
    fprintf('\n');
end
end


function [w] = ItotIndicator(d,Icritical)
w = cell(d,1);
w{1} = zeros(1,3,d-Icritical+1);
w{1}(1,:,1) = [0 1 0];
w{1}(1,:,2) = [1 0 1];
for i=2:d-1
    w{i} = zeros(d-Icritical+1, 3, d-Icritical+1);
    for j=1:(d-Icritical)
        w{i}(j,:,j) = [0 1 0];
        w{i}(j,:,j+1) = [1 0 1];
    end
    w{i}(d-Icritical+1,:,d-Icritical+1) = [0 1 0];
end
w{d} = ones(d-Icritical+1, 3);
w{d}(d-Icritical+1, :) = [0 1 0];
w = cell2core(tt_tensor, w);
end


