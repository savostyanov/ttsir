% Simulate new data for network inference or load it from existing file

% SSA - compute exact samples for the given network - these will be observed
Nsim =  parse_parameter('Number of trajectories to observe (0 to use previous data)', 0);
if (Nsim>0)
    % Populate z and w from adjacency matrix
    % z is d x M
    z = zeros(d, 2*d+numel(iadj));
    wfun = cell(2, 2*d+numel(iadj));
    for m=1:numel(iadj)
        % infection iadj -> jadj
        z(jadj(m),m) = 1; % increment the state by 1: only 0->1 transition is possible, exactly what we need
        wfun{1,m} = @(x)beta(min(m,numel(beta))) * double(x(:,1)==0).*double(x(:,2)==1);
        wfun{2,m} = [jadj(m), iadj(m)];
    end
    % SI recovery
    for m=1:d
        z(m,numel(iadj)+m) = -1;
        wfun{1,numel(iadj)+m} = @(x)gamma(min(m,numel(gamma))) * double(x==1);
        wfun{2,numel(iadj)+m} = m;
    end
    % Bath infection
    for m=1:d
        z(m,d+numel(iadj)+m) = 1;
        wfun{1,d+numel(iadj)+m} = @(x)delta(min(m,numel(delta))) * double(x==0);
        wfun{2,d+numel(iadj)+m} = m;
    end

    % Run this for SSA
    ttimes_ssa = tic;
    Xobs = ssa(z,wfun,x0,T,Nsim,tfix, 100);
    ttimes_ssa = toc(ttimes_ssa);
    ind0 = randi(2, 1, d*(d-1)/2);
    ind1 = scoring_initial_net(Xobs);
    err0 = norm(ind0 - imax_ex,1);
    nnz0 = sum(ind0==2);
    err1 = norm(ind1 - imax_ex,1);
    nnz1 = sum(ind1==2);
    save(sprintf('Xobs-Net%d-d%d-irun%d.mat', network_type, d, irun), 'Xobs', 'tfix', 'ind0', 'd', 'imax_ex');
    figure(1);
    scatter3(reshape(repmat(tfix',d,1),[],1), reshape(repmat((1:d)',1,numel(tfix)),[],1), reshape(mean(Xobs,1),[],1), 20, reshape(mean(Xobs,1),[],1), 'filled'); view(2);
    colorbar;
    fprintf('\terr\tnnz\n0\t%d\t%d\n1\t%d\t%d\n', err0, nnz0, err1, nnz1);
else
    load(sprintf('Xobs-Net%d-d%d-irun%d.mat', network_type, d, irun));
end
ind0 = scoring_initial_net(Xobs);

L_ex = cme_likelihood_si_bath(imax_ex, beta, gamma, delta, tfix, tol, Xobs)
