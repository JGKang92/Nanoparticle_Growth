function trange = Determine_trange(rc, kappa_c, lambda, sigmac)
%% Parameters
% x(1): sigma_c
% x(2): kappa_c
% x(3): a1s (variance of S(r))
% x(4): l (ligand length)
rmax = 300;
tmin = 8;
tmax = 180;
xlen = 5000;
lent = 5000;
a1s = (9.93*1e-3).^2;   % v^2 value in manuscript
l = 2.87;               % Passivative ligand length
b1k = 7.64;             % gamma value in manuscript
%% Sigma related function
sigmat = @(t) lambda.*(1.7773331025798687 + 0.04422044832970424.*t - 0.0008212384775235609.*t.*t + ...
    7.651821072767395.*(10.^-6).*t.*t.*t - 2.153266857581242.*(10.^-8).*t.*t.*t.*t);
dsigmat = @(t) lambda.*(0.04422044832970424 - 0.0016424769550471218.*t + 0.000022955463218302186.*t.*t - ...
    8.613067430324968.*(10.^-8).*t.*t.*t);
%% Kappa and Sink function
% c1k = 3.007832240727688;
c1k = sigmac;
kappa_t = @(kappa_c,t) kappa_c./(exp(b1k.*(sigmat(t) - c1k))+1);
Sr = @(r,t) exp(-(r-(sigmat(t)+l)).^2/(2*a1s));
%% Initial distribution
Pinit = @(sigmac,r) (r >= sigmac & r <= 270);
%% Relative diffusion coefficient
% This file is given in Dr.m
Drt = @(t) Dr(t, lambda);
%% Solve PDE
pdefun = @(x,t,u,dudx) coalescence_pde(x, t, u, dudx, Drt, sigmat, dsigmat, kappa_t, Sr, kappa_c);
icfun = @(x) coalescence_initial_cond(x, Pinit, sigmat, rc);
bcfun = @(xl,ul,xr,ur,t) coalescence_bd_cond(xl,ul,xr,ur,t,sigmat);
xmesh = linspace(0,rmax,xlen);
tspan = linspace(tmin,tmax,lent);
sol = pdepe(0,pdefun,icfun,bcfun,xmesh,tspan);
[colnum, ~] = size(sol);
if colnum ~= lent
    error('Column number is not satisfied.')
end
%% kappa_f
kappa_fval = zeros(1,lent);
for idx=1:lent
    tval = tspan(idx);
    ufun_pp = csapi(xmesh, sol(idx,:));
    ufun = @(x) 4*pi*(x+sigmat(tval)).^2.*kappa_t(kappa_c,tval).*Sr(x+sigmat(tval),sigmat(tval)).*ppval(ufun_pp,x);
    kappa_fval(idx) = integral(ufun,xmesh(1),xmesh(end));
end
kappa_f_pp = csapi(tspan, kappa_fval);
kappa_f = @(t) ppval(kappa_f_pp,t);
%% Survival probability
survt = zeros(1,length(tspan));
for idx=1:length(tspan)
    survt(idx) = 1./(1 + integral(@(t)kappa_f(t),tspan(1),tspan(idx)));
end
%% Coalescence time distribution
wtd_theory = survt.^2.*kappa_f(tspan);
wtd_pp = csapi(tspan,wtd_theory);
wtd_fun = @(x)ppval(wtd_pp,x);
norm_const = integral(wtd_fun,tspan(1),tspan(end));
wtd_fun = @(x)ppval(wtd_pp,x) ./ norm_const;
%% Find interquartile range for coalescence time distribution
q1_bool = true;
trange = zeros(1,2);
for idx=1:lent
    intval = integral(wtd_fun, tspan(1), tspan(idx));
    if (intval >= 0.25 && q1_bool)
        q1_bool = false;
        trange(1) = tspan(idx);
    elseif (intval >= 0.75)
        trange(2) = tspan(idx);
        break;
    end
end
if (trange(2) == 0)
    warning('Zero value')
end