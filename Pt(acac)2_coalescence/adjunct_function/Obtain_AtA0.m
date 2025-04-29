function [tval, AtA0_theory] = Obtain_AtA0(x, AtA0_exp, str)
%% Parameters
% x(1): sigma_c
% x(2): kappa_c
% x(3): a1s (variance of S(r))
% x(4): l (ligand length)
% x(5): P0 value
rmax = 200;
tmax = 180;
sigmac = x(1);
kappa_c = x(2);
a1s = x(3);
l = x(4);
xlen = 5000;
lent = 5000;
tmin = 8;
% P0 = 0.00394;
% P0tval = 93;
P0 = x(5);
P0tval = x(6);
b1k = x(7);
%% Sigma related function
sigmat = @(t) 1.7773331025798687 + 0.04422044832970424.*t - 0.0008212384775235609.*t.*t + ...
    7.651821072767395.*(10.^-6).*t.*t.*t - 2.153266857581242.*(10.^-8).*t.*t.*t.*t;
dsigmat = @(t) 0.04422044832970424 - 0.0016424769550471218.*t + 0.000022955463218302186.*t.*t - ...
    8.613067430324968.*(10.^-8).*t.*t.*t;
%% Kappa and Sink function
% b1k = 0.0683700205567837;
c1k = 3.007832240727688;
if strcmp(str, 'cdf')
%     kappa_t = @(kappa_c, t) kappa_c.*(0.5.*(2-erf(c1k./sqrt(2.*b1k))+erf((-sigmat(t)+c1k)./sqrt(2.*b1k))));
%     kappa_t = @(kappa_c,t) 0.69545;
    kappa_t = @(kappa_c,t) kappa_c./(exp(b1k.*(sigmat(t) - c1k))+1);
%     kappa_t = @(kappa_c,t) kappa_c;
elseif strcmp(str, 'pdf')
    kappa_t = @(kappa_c, t) kappa_c.*exp(-(sigmat(t)-c1k).^2/(2.*b1k))./sqrt(2.*pi.*b1k);
else
    error('str should be one of cdf or pdf');
end
% Sr = @(r,t) exp(-(r-(sigmat(t)+l)).^2/(2*a1s))./(2.*pi.*a1s);
Sr = @(r,t) exp(-(r-(sigmat(t)+l)).^2/(2*a1s));
% savedata = [linspace(0,10,100)', kappa_c./(exp(b1k.*(linspace(0,10,100)' - c1k))+1)];
% save('kappa_t_value.dat','savedata','-ascii')
%% Initial distribution
Pinit = @(sigmac,r) (r >= sigmac & r <= 270);
% Pinit = @(sigmac,r) (r>=1.7773331025798687 & r<=270);
%% Relative diffusion coefficient
% This file is given in Dr.m
%% Solve PDE
pdefun = @(x,t,u,dudx) coalescence_pde(x, t, u, dudx, @Dr, sigmat, dsigmat, kappa_t, Sr, kappa_c);
icfun = @(x) coalescence_initial_cond(x, Pinit, sigmat, sigmac);
bcfun = @(xl,ul,xr,ur,t) coalescence_bd_cond(xl,ul,xr,ur,t,sigmat);
xmesh = linspace(0,rmax,xlen);
tspan = linspace(tmin,tmax,lent);
sol = pdepe(0,pdefun,icfun,bcfun,xmesh,tspan);
pde_data = [xmesh', sol'];
t_data = tspan';
save('pair_corrleation_dist.dat','pde_data','-ascii');
save('pair_correlation_time.dat','t_data','-ascii');
%% kappa_f value
kappa_fval = zeros(1,lent);
for idx=1:lent
    tval = tspan(idx);
    ufun_pp = csapi(xmesh, sol(idx,:));
    ufun = @(x) 4*pi*(x+sigmat(tval)).^2.*kappa_t(kappa_c,tval).*Sr(x+sigmat(tval),sigmat(tval)).*ppval(ufun_pp,x);
    kappa_fval(idx) = integral(ufun,xmesh(1),xmesh(end));
end
kappa_f_pp = csapi(tspan, kappa_fval);
kappa_f = @(t) ppval(kappa_f_pp,t);
figure(2)
clf
plot(tspan, kappa_f(tspan), '-')
savedata = [tspan', kappa_f(tspan)'];
save('kf_theory_umax_200.dat','savedata','-ascii');
%% [A](t)/[A](0) value
% Solve the ODE for bimolecular reaction with time-dependent rate
% coefficient, dSa/dt = p0 - K_f(t)Sa^2
kinetic_eq = @(t,y) P0.*(t< P0tval) - kappa_f(t).*y.*y;
yinit = AtA0_exp(1);
[tval, AtA0_theory] = ode45(kinetic_eq,tspan, yinit);