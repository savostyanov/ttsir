% Plot statistics of the network configuration error
check_tt;
addpath('inferdata');

d = parse_parameter('dimension (number of nodes)', 9);
network_type = parse_parameter('Network type (11 - chain, 14 - Austria)', 11);

experiments = 1:42; % Which experiments to include

% Mean+-std of index vs. evals
figure(4);
evalgrid = [1:4 5:5:400];  % (Nepochs-1)*d*(d-1)/2
clear ierr ierr2
for irun=experiments
    load(sprintf('Net%d-d%d-irun%d-methods%d.mat', network_type, d, irun, 2));
    ierr2(irun,:) = interp1(eval_mcmc2(1:end-1), sum(abs(imax_mcmc2-imax_ex),2), evalgrid, 'previous', 'extrap');
end
for irun=experiments
    load(sprintf('Net%d-d%d-irun%d-methods%d.mat', network_type, d, irun, 4));
    ierr(irun,:) = interp1(eval_mcmc(1:end-1), sum(abs(imax_mcmc-imax_ex),2), evalgrid, 'previous', 'extrap');
end
ierr = ierr(experiments, :);
ierr2 = ierr2(experiments, :);
mierr = mean(ierr, 'omitnan');
lierr = mierr - std(ierr, 'omitnan'); % quantile(ierr, 0.1);    
hierr = mierr + std(ierr, 'omitnan'); % quantile(ierr, 0.9);
mierr2 = mean(ierr2, 'omitnan');
lierr2 = mierr2 - std(ierr2, 'omitnan'); % quantile(ierr2, 0.1);
hierr2 = mierr2 + std(ierr2, 'omitnan'); % quantile(ierr2, 0.9);

hold off;
fill([evalgrid fliplr(evalgrid)], [hierr fliplr(lierr)], 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold on
plot(evalgrid, mierr, 'g', 'LineWidth', 2);
fill([evalgrid fliplr(evalgrid)], [hierr2 fliplr(lierr2)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
plot(evalgrid, mierr2, 'r', 'LineWidth', 2);

hold off;
legend('mcmc mean+-std', 'mcmc mean', 'mcmc2 mean+-std', 'mcmc2 mean')
xlabel('N evaluations')
title('|i-i_*|_1')

dat = [evalgrid', mierr'];
save(sprintf('Net%d-d%d-grid-mcmc-mean.dat', network_type, d), '-ascii', 'dat');
dat = [[evalgrid fliplr(evalgrid)]', [hierr fliplr(lierr)]'];
save(sprintf('Net%d-d%d-grid-mcmc-pmsd.dat', network_type, d), '-ascii', 'dat');

dat = [evalgrid', mierr2'];
save(sprintf('Net%d-d%d-grid-mcmc2-mean.dat', network_type, d), '-ascii', 'dat');
dat = [[evalgrid fliplr(evalgrid)]', [hierr2 fliplr(lierr2)]'];
save(sprintf('Net%d-d%d-grid-mcmc2-pmsd.dat', network_type, d), '-ascii', 'dat');
