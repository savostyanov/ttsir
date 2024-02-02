% MCMC algorithms for network inference

% MCMC pre-shuffled 
if (bitand(methods, 2) == 2)
    touter = tic;
    L_mcmc2 = zeros(Nepochs*d*(d-1)/2, 1);
    ind_mcmc2 = zeros(numel(L_mcmc2), d*(d-1)/2);
    imax_mcmc2 = ind0;
    ind_mcmc2(1,:) = ind0;
    Lmax_mcmc2 = cme_likelihood_si_bath(ind_mcmc2(1,:), beta, gamma, delta, tfix, tol, Xobs)
    L_mcmc2(1) = Lmax_mcmc2;
    ttimes_mcmc2 = 0;
    eval_mcmc2 = 1;
    Nrej2 = 0;
    i = 2;
    for m=1:Nepochs
        [~,i_shuffle] = sort(rand(1,d*(d-1)/2))
        for j=1:d*(d-1)/2
            ind_new = ind_mcmc2(i-1,:);
            ind_new(i_shuffle(j)) = 3 - ind_new(i_shuffle(j));
            L_new = cme_likelihood_si_bath(ind_new, beta, gamma, delta, tfix, tol, Xobs);
            alpha = L_new - L_mcmc2(i-1);
            alpha = exp(alpha);                    % No tempering
            if (alpha<rand)
                % reject, continue with (i-1)th data
                ind_mcmc2(i,:) = ind_mcmc2(i-1,:);
                L_mcmc2(i) = L_mcmc2(i-1);
                Nrej2 = Nrej2 + 1;
            else
                ind_mcmc2(i,:) = ind_new;
                L_mcmc2(i) = L_new;
                if (L_new > Lmax_mcmc2(end))
                    imax_mcmc2(end+1,:) = ind_new
                    Lmax_mcmc2(end+1) = L_new
                    ttimes_mcmc2(end+1) = toc(touter);
                    eval_mcmc2(end+1) = i
                    save(sprintf('Net%d-d%d-irun%d-methods%d.mat', network_type, d, irun, methods));
                end
            end
            i = i+1;
        end
    end
    fprintf('MCMC2 completed with %g%% rejections\n', Nrej2/numel(L_mcmc2)*100);
    
    W_mcmc2 = full(ind_to_adj(d, imax_mcmc2(end,:)))
    ttimes_mcmc2(end+1) = toc(touter)
    eval_mcmc2(end+1) = i;
    
    err_W_mcmc2 = W_mcmc2 - W_ex
    
    save(sprintf('Net%d-d%d-irun%d-methods%d.mat', network_type, d, irun, methods));    
end

% MCMC optimizer - proposal is flipping one uniformly random index at a time
if (bitand(methods, 4) == 4)
    touter = tic;
    L_mcmc = zeros(Nepochs*d*(d-1)/2, 1);
    ind_mcmc = zeros(numel(L_mcmc), d*(d-1)/2);
    imax_mcmc = ind0;
    ind_mcmc(1,:) = ind0;
    Lmax_mcmc = cme_likelihood_si_bath(ind_mcmc(1,:), beta, gamma, delta, tfix, tol, Xobs)
    L_mcmc(1) = Lmax_mcmc;
    ttimes_mcmc = 0;
    eval_mcmc = 1;
    Nrej = 0;
    for i=2:numel(L_mcmc)
        ind_new = ind_mcmc(i-1,:);
        i_flip = randi(d*(d-1)/2);
        ind_new(i_flip) = 3 - ind_new(i_flip);
        L_new = cme_likelihood_si_bath(ind_new, beta, gamma, delta, tfix, tol, Xobs);
        alpha = L_new - L_mcmc(i-1);
        alpha = exp(alpha);                    % No tempering
        if (alpha<rand)
            % reject, continue with (i-1)th data
            ind_mcmc(i,:) = ind_mcmc(i-1,:);
            L_mcmc(i) = L_mcmc(i-1);
            Nrej = Nrej + 1;
        else
            ind_mcmc(i,:) = ind_new;
            L_mcmc(i) = L_new;
            if (L_new > Lmax_mcmc(end))
                imax_mcmc(end+1,:) = ind_new
                Lmax_mcmc(end+1) = L_new
                ttimes_mcmc(end+1) = toc(touter);
                eval_mcmc(end+1) = i
                save(sprintf('Net%d-d%d-irun%d-methods%d.mat', network_type, d, irun, methods));
            end
        end
    end
    fprintf('MCMC completed with %g%% rejections\n', Nrej/numel(L_mcmc)*100);
    
    W_mcmc = full(ind_to_adj(d, imax_mcmc(end,:)))
    ttimes_mcmc(end+1) = toc(touter)
    eval_mcmc(end+1) = i;
    
    err_W_mcmc = W_mcmc - W_ex
    
    save(sprintf('Net%d-d%d-irun%d-methods%d.mat', network_type, d, irun, methods));
end
