% Plot statistics of the log-likelihood error
check_tt;
addpath('inferdata');

d = parse_parameter('dimension (number of nodes)', 9);
network_type = parse_parameter('Network type (11 - chain, 14 - Austria)', 11);

experiments = 1:42; % Which experiments to include

% Mean+-std of log-likelihood vs. evals
figure(3);
evalgrid = [1:4 5:5:400];
clear Ldiff Ldiff2
for irun=experiments
    load(sprintf('Xobs-Net%d-d%d-irun%d.mat', network_type, d, irun));
    load(sprintf('Net%d-d%d-irun%d-methods%d.mat', network_type, d, irun, 2));
    Ldiff2(irun,:) = interp1(eval_mcmc2(1:end-1), Lmax_mcmc2, evalgrid, 'previous', 'extrap') - L_ex;
end
for irun=experiments
    load(sprintf('Xobs-Net%d-d%d-irun%d.mat', network_type, d, irun));
    load(sprintf('Net%d-d%d-irun%d-methods%d.mat', network_type, d, irun, 4));
    Ldiff(irun,:) = interp1(eval_mcmc(1:end-1), Lmax_mcmc, evalgrid, 'previous', 'extrap') - L_ex;
end
Ldiff = Ldiff(experiments, :);
Ldiff2 = Ldiff2(experiments, :);
mLdiff =      mean(Ldiff, 'omitnan');
lLdiff =  mLdiff - std(Ldiff, 'omitnan'); % quantile(Ldiff, 0.1);
hLdiff =  mLdiff + std(Ldiff, 'omitnan'); % quantile(Ldiff, 0.9);
mLdiff2 =     mean(Ldiff2, 'omitnan');
lLdiff2 = mLdiff2 - std(Ldiff2, 'omitnan'); % quantile(Ldiff2, 0.1);
hLdiff2 = mLdiff2 + std(Ldiff2, 'omitnan'); % quantile(Ldiff2, 0.9);
hold off;

fill([evalgrid fliplr(evalgrid)], [hLdiff fliplr(lLdiff)], 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
hold on
plot(evalgrid, mLdiff, 'g', 'LineWidth', 2);
fill([evalgrid fliplr(evalgrid)], [hLdiff2 fliplr(lLdiff2)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
plot(evalgrid, mLdiff2, 'r', 'LineWidth', 2);

hold off;
legend('mcmc mean+-std', 'mcmc mean', 'mcmc2 mean+-std', 'mcmc2 mean')
xlabel('N evaluations')
title('L-L_*')


dat = [evalgrid', mLdiff'/log(10)];
save(sprintf('Net%d-d%d-like-mcmc-mean.dat', network_type, d), '-ascii', 'dat');
dat = [[evalgrid fliplr(evalgrid)]', [hLdiff fliplr(lLdiff)]'/log(10)];
save(sprintf('Net%d-d%d-like-mcmc-pmsd.dat', network_type, d), '-ascii', 'dat');

dat = [evalgrid', mLdiff2'/log(10)];
save(sprintf('Net%d-d%d-like-mcmc2-mean.dat', network_type, d), '-ascii', 'dat');
dat = [[evalgrid fliplr(evalgrid)]', [hLdiff2 fliplr(lLdiff2)]'/log(10)];
save(sprintf('Net%d-d%d-like-mcmc2-pmsd.dat', network_type, d), '-ascii', 'dat');
