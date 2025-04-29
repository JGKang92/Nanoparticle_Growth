% Partial differential equation solver (20210201 version)
clear
addpath('functions')
%% Import experimental data
AtA0_exp = load('AtA0value.dat');
%% Perform fitting
% x(1): r_c
% x(2): kappa_c
% x(3): a1s (variance of S(r), nu)
% x(4): l (ligand length)
% P0 = x(5) Initial slope
% P0tval = x(6) cutoff value (gamma)
% x0 = [44,11.6229739205327,0.00960299780366661,2.84967169163257,0.00367433442774329,96,7.60142298460863];
x0 = [44, 11.6, 0.00993^2, 2.87, 0.1, 96, 7.64];
fitfun = @(x)PDE_based_costfun(x, AtA0_exp, 100);
lb = x0*0.5;
ub = x0*1.5;
lb(4) = 2.7;
ub(4) = 4;
nonlcon = @(x) PDE_based_constraint(x, AtA0_exp);
% options = optimoptions('fmincon', 'Display', 'iter', 'UseParallel', true, 'StepTolerance', 1e-22, 'ConstraintTolerance', 1e-16);
options = optimoptions('fmincon', 'Display', 'iter', 'UseParallel', true);
xfit = fmincon(fitfun, x0, [], [], [], [], lb, ub, nonlcon, options);
%% Plot the result
[tval, AtA0_theory] = Obtain_AtA0(xfit, AtA0_exp);
figure(1)
clf
plot(8:180, AtA0_exp, '.', tval, AtA0_theory, '-')
%%
% savedata = [(8:180)', AtA0_exp'];
% save('AtA0_value_exp.dat','savedata','-ascii')
%% Draw figure
load('coal_time.mat')
[N, edges] = histcounts(coal_time, 5, 'Normalization','pdf');
figure(3)
clf
subplot(1,2,1)
histogram(coal_time,5,'Normalization','pdf')
hold on
plot(tval, wtd_theory./norm_const, '-', 'MarkerSize', 20)
hold off
subplot(1,2,2)
plot(1:180, tot_radi_mean, '.', 8:180, movmean(tot_mean_theory,10), '-')
savedata = [tval, (wtd_theory./norm_const)'];
save('coalescence_time_dist_theory.dat', 'savedata', '-ascii');
%% Experimental kf(t)
P0 = xfit(5);
P0tval = xfit(6);
tot_radi_mean_pp = csapi(8:180, AtA0_exp);
tot_radi_mean_f = @(t) ppval(tot_radi_mean_pp,t);
d_tot_radi_mean_pp = fnder(tot_radi_mean_pp,1);
d_tot_radi_mean_f = @(t) fnval(d_tot_radi_mean_pp,t);
kf_fun = @(t) (P0.*(t< P0tval) - d_tot_radi_mean_f(t))./tot_radi_mean_f(t).^2;
data = [(8:180)', movmean(kf_fun(8:180),20)'];
figure(4)
clf
plot(8:180, movmean(kf_fun(8:180),20))
save('k_f_exp.dat','data','-ascii');