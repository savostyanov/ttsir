%   [outputs] = EXPECTATIONS_SSA(xfix, outfun)
% Compute average values of functions of interest evaluated at SSA samples
% Inputs:
%   xfix: Nsim x d x numel(tfix) array of copy numbers returned from SSA
%   outfun: l x R cell array of QoI functions (handles), l=1,2.
%           If l==1: each function outfun{i} takes a I x d double array of 
%           X values and returns a I x 1 double array of the function values, 
%           where I can be any natural number. In other words, X(j,:) is a
%           vector of copy numbers of all species (the state), where the 
%           propensity function should be evaluated.
%           If l==2: outfun{2,i} contains a vector of coordinates the
%           function outfun{1,i} depends on. In this case, outfun{1,i} should
%           take a I x numel(outfun{2,i}) array of X{outfun{2,i}} values, and
%           return a I x 1 vector of the values.
%
%           You can use the function restore_default_outfun to produce the
%           outfun array for computing the first and second moments.
% Outputs:
%   outputs: numel(tfix) x R array of expectation values

function [outputs] = expectations_ssa(xfix, outfun)
[Nsim,d,nfix] = size(xfix);

% Determine the particular format of outfun
if (isa(outfun,'cell'))
    if (size(outfun,1)==2)
        dmap = outfun(2,:);
        outfun = outfun(1,:);
    else
        dmap = num2cell(repmat(1:d, numel(outfun), 1), 2);
    end
end

outputs = zeros(nfix, numel(outfun));
% Loop over time steps
for it=1:nfix
    if (isa(outfun, 'cell'))
        out_local = zeros(Nsim, numel(outfun));
        for i=1:numel(outfun)
            out_local(:,i) = feval(outfun{i}, xfix(:, dmap{i}, it));
        end
    else
        out_local = feval(outfun, xfix(:, :, it));
    end
    outputs(it, :) = mean(out_local,1);
end
end
