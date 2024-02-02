function [L] = cme_likelihood_si_bath(ind, beta, gamma, delta, tfix, tol, Xobs)
% ind \in {1,2}^{d(d-1)/2} parametrises off-diag entries of W
% Xobs(isample, idimension, itime) \in {0,1}
d = size(Xobs,2);
M = size(ind,1);
assert(size(ind,2)==d*(d-1)/2, 'ind should contain d*(d-1)/2 parameters')
assert(size(Xobs,3)==numel(tfix), '3d dim of Xobs should be the number of time points in tfix')

method = 0; % 0: TT, 1: SSA

nt = 12; % Chebyshev degree
Nsim = 1000; % N SSA trajectories tuneable

if (method==0)
    % Matrix parts
    Afix = cell(1,d);
    % Recovery and bath reinfection
    for i=1:d
        Ai = ([0,1;1,0] - eye(2)) * diag([0;1]) * gamma + ([0,1;1,0] - eye(2)) * diag([1;0]) * delta;
        Ai = tt_matrix(Ai);
        if (i>1)
            Ai = tkron(tt_eye(2,i-1), Ai);
        end
        if (i<d)
            Ai = tkron(Ai, tt_eye(2,d-i));
        end
        Afix{1,i} = tt_tensor(Ai);
    end
    Afix = amen_sum(Afix, ones(numel(Afix),1), 1e-12, 'y0', tt_tensor(tt_eye(2,d)), 'verb', 0);
    Afix = tt_matrix(Afix);
end

L = zeros(M,1); % storage for log-likelihoods
% Loop over cross samples
for m=1:M
    W = ind_to_adj(d, ind(m,:));
    % Try permutations
    [VV,LL] = eigs(diag(sum(W)) - W, 2, 'SA');
    [~,ii] = max(diag(LL));
    [~,p] = sort(VV(:,ii));
    W = W(p,p);
    Xobs = Xobs(:, p, :);
    
    [iadj, jadj, ~] = find(W);

    if (method==1)  % SSA

        % Populate z and w from adjacency matrix
        % z is d x M
        z = zeros(d, 2*d+numel(iadj));
        wfun = cell(2, 2*d+numel(iadj));
        for r=1:numel(iadj)
            % infection iadj -> jadj
            z(jadj(r),r) = 1; % increment the state by 1: only 0->1 transition is possible, exactly what we need
            wfun{1,r} = @(x)beta * double(x(:,1)==0).*double(x(:,2)==1);
            wfun{2,r} = [jadj(r), iadj(r)];
        end
        % recovery
        for r=1:d
            z(r,numel(iadj)+r) = -1; % decrement the state by 1: only 1->0 transition is possible, exactly what we need
            wfun{1,numel(iadj)+r} = @(x)gamma * double(x==1);
            wfun{2,numel(iadj)+r} = r;
        end
        % Bath infection
        for r=1:d
            z(r,d+numel(iadj)+r) = 1; 
            wfun{1,d+numel(iadj)+r} = @(x)delta * double(x==0);
            wfun{2,d+numel(iadj)+r} = r;
        end
        

    else            % TT

        % Infection
        Acell = cell(1,numel(iadj));
        for i=1:numel(iadj)
            % iadj -> jadj
            Ai = beta*diag([0;1]);
            Aj = ([0,1;1,0] - eye(2)) * diag([1;0]);
            Ai = tt_matrix(Ai);
            Aj = tt_matrix(Aj);
            if (iadj(i)<jadj(i))
                if (iadj(i)>1)
                    Ai = tkron(tt_eye(2,iadj(i)-1), Ai);
                end
                if (jadj(i)<d)
                    Aj = tkron(Aj, tt_eye(2, d-jadj(i)));
                end
                if (jadj(i)-iadj(i)>1)
                    Ai = tkron(Ai, tt_eye(2, jadj(i)-iadj(i)-1));
                end
                Ai = tkron(Ai, Aj);
            else
                if (jadj(i)>1)
                    Aj = tkron(tt_eye(2,jadj(i)-1), Aj);
                end
                if (iadj(i)<d)
                    Ai = tkron(Ai, tt_eye(2, d-iadj(i)));
                end
                if (iadj(i)-jadj(i)>1)
                    Aj = tkron(Aj, tt_eye(2, iadj(i)-jadj(i)-1));
                end
                Ai = tkron(Aj, Ai);
            end
            Acell{1,i} = tt_tensor(Ai);
        end

        A = Afix;
        if (~isempty(iadj))
            Ainf = amen_sum(Acell, ones(numel(iadj),1), 1e-12, 'y0', tt_tensor(tt_eye(2,d)), 'verb', 0);
            Ainf = tt_matrix(Ainf);
            A = round(Afix + Ainf, 1e-12);
        end
        A = core2cell(A);

    end

    Lloc = nan(numel(tfix)-1, size(Xobs,1));
    ranks = nan(numel(tfix)-1, size(Xobs,1));
    ttt = tic;
    % Find SSA trajectories and probability over them
    for j=1:size(Xobs,1)
        for i=2:numel(tfix)
            if (method==1)  % SSA
                Xrand = ssa(z,wfun,Xobs(j,:,i-1),tfix(i)-tfix(i-1),Nsim,tfix(i)-tfix(i-1), Nsim+1); % Xrand of size Nsim x d
                Lloc(i-1,j) = log(sum(double(all(Xrand==Xobs(j,:,i), 2)))/Nsim);
            else            % TT
                u = arrayfun(@(i)double((0:1)==i), Xobs(j,:,i-1)', 'uni', 0);
                Uglobal = [u; {ones(1,nt)}];
                A{d} = A{d} * (tfix(i)-tfix(i-1));
                [Uglobal,tglobal] = tamen(Uglobal,A,tol,struct('time_error_damp', 0, 'verb', 0));
                A{d} = A{d} / (tfix(i)-tfix(i-1));
                ranks(i-1,j) = max(cellfun(@(u)size(u,1), Uglobal));
                u = extract_snapshot(Uglobal, tglobal, 1);
                u = cellfun(@(ui)reshape(ui, size(ui,1), size(ui,2), size(ui,4)), u, 'uni', 0);
                Li = tt_sample_ind(u, Xobs(j,:,i)+1);
                Lloc(i-1,j) = sum(log(abs(Li)));
            end
        end
    end

    Lloc = Lloc(:);
    Lloc(isnan(Lloc)) = [];
    L(m) = sum(Lloc);

    fprintf('LValue %d/%d took %g sec.  %g\n', m, M, toc(ttt), L(m));
end
end
