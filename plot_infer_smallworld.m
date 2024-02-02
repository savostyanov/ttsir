% Plot log-likelihood errors in small world experiments

d = 15;
network_type = 17;
experiments = 1:42;

% Plot individual 2D loglikelihoods
figure(1);
err_L = zeros(d,d);
for irun=experiments
    load(sprintf('Net%d-d%d-irun%d.mat', network_type, d, irun));
    load(sprintf('Xobs-Net%d-d%d-irun%d.mat', network_type, d, irun));
    err_L = err_L + L - L(1,8);
    figure(1);
    imagesc(L);
    [~,imax] = max(L(:));
    [imax,jmax] = ind2sub([d d], imax);
    title(sprintf('log-likelihood, %dth run, MLE=[%d %d]', irun, imax, jmax));
    xlabel('i');
    ylabel('j');
    colorbar;
    drawnow;
end
err_L = err_L / numel(experiments);
figure(2);
imagesc(err_L);
title('E[L-L_*]');
xlabel('i');
ylabel('j');
colorbar;

