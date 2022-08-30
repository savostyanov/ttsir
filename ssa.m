%   [xfix] = SSA(z,wfun,x0,T, Nsim, tfix, print_tick)
% Vanilla Gillespie's Stochastic Simulation Algorithm.
% Inputs:
%     z: d x R stoichiometry matrix, where d is the number of species, and
%           R is the number of reactions.
%     wfun: l x R cell array of propensity functions (handles), l=1,2.
%           If l==1: each function wfun{i} takes a I x d double array of 
%           X values and returns a I x 1 double array of W values, 
%           where I can be any natural number. In other words, X(j,:) is a
%           vector of copy numbers of all species (the state), where the 
%           propensity function should be evaluated.
%           If l==2: wfun{2,i} contains a vector of coordinates the
%           function wfun{1,i} depends on. In this case, wfun{1,i} should
%           take a I x numel(wfun{2,i}) array of X{wfun{2,i}} values, and
%           return a I x 1 vector of W values.
%     x0: a d x 1 double vector of deterministic initial state values.
%     T: time interval
%     Nsim: number of simulated trajectories
%     tfix: an array of time points where the solution should be evaluated
%     print_tick: (optional) A step (in %) in the total progress after 
%       which the progress is printed. Default value is 10.
% Outputs:
%   xfix: Nsim x d x numel(tfix) array of copy numbers.
function [xfix]=ssa(z,wfun,x0, T, Nsim, tfix, print_tick, outfun)

[d,R] = size(z);
nfix = numel(tfix);
tfix = reshape(tfix, 1, nfix);
x0 = reshape(x0, 1, d);
if (nargin<8)||(isempty(outfun))
    xfix = zeros(Nsim, d, nfix);
    Rout = 0;
else
    Rout = size(outfun,2);
    if (size(outfun,1)==2)
        dmapo = outfun(2,:);
        outfun = outfun(1,:);
    else
        dmapo = num2cell(repmat(1:d, Rout, 1), 2);
    end
    xfix = zeros(Rout, nfix);
end

if (nargin<7)
    print_tick = 100;
end

% Determine the particular format of wfun
if (size(wfun,1)==2)
    dmapw = wfun(2,:);
    wfun = wfun(1,:);
else
    dmapw = num2cell(repmat(1:d, R, 1), 2);
end

% Loop over trajectories (inefficient but the only working way now)
for isim=1:Nsim
    xhist = zeros(1000,d);
    thist = zeros(1000,1);
    x = x0;
    xhist(1,:) = x;
    t = 0;    
    tstep = 1;  
    % Run SSA trajectory
    while (t<=T)
        % compute propensities at the current state
        w = zeros(1, R);
        for i=1:R
            w(:,i) = feval(wfun{i}, x(:,dmapw{i}));
        end
        w_recp = sum(w, 2);
        w_recp = 1./w_recp;
        % new time
        r = rand(1,1);
        tau = log(1./r).*w_recp;
        t = t + tau;
        % new reaction
        w_recp = repmat(w_recp, 1, R);
        w = w.*w_recp;
        w = cumsum(w,2); % ranges from 0 to 1
        r = rand(1, 1);
        react_indices = sum(r>=w, 2)+1; % vectorised selection
        % update x here
        x = x + z(:,react_indices)';
        % Save x to the history
        tstep = tstep + 1;
        if (numel(thist)<tstep)
            % Allocate more memory if necessary
            thist = [thist; zeros(1000,1)];
            xhist = [xhist; zeros(1000,d)];            
        end
        thist(tstep) = t;
        xhist(tstep,:) = x;
    end
    thist = thist(1:tstep);
    xhist = xhist(1:tstep, :);
    
    % Parse the trajectory into time steps
    if (Rout==0)
        for j=1:nfix
            tstep = find(thist>tfix(j), 1)-1;
            xfix(isim,:,j) = xhist(tstep, :);
        end
    else
        xfixlocal = zeros(Rout, nfix);
        for j=1:nfix
            tstep = find(thist>tfix(j), 1)-1;
            xt = xhist(tstep, :);
            for i=1:Rout
                xfixlocal(i,j) = feval(outfun{i}, xt(dmapo{i}));
            end
        end
        xfix = xfix + xfixlocal;
    end
    if (mod(isim, print_tick)==0)
        fprintf('SSA trajectory %d of %d done\n', isim, Nsim);
    end
end

if (Rout>0)
    xfix = xfix.'/Nsim; % Return nfix x R for compatibility with tt_cme
end
end
