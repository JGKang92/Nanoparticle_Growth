function [c,f,s] = coalescence_pde(x,t,u,dudx,Dr, sigmat, dsigmat, kappa_t,Sr, kappa_c)
    % Dr: relative diffusion coefficient, in function handle
    % sigmat: sigma value in function handle
    % dsigmat: time derivative of simgat value in function handle
    % kappa_t: time-dependent reaction rate coefficient in function handle
    % Sr: distance-dependent sink function in function handle
    % Note: 3-dimensional case, PDE with pair-correlation!
    c = 1./Dr(t);
    f = dudx;
    s = (2./(x+sigmat(t)) + dsigmat(t)./Dr(t))*dudx - kappa_t(kappa_c, t).*Sr(x + sigmat(t),t).*u./Dr(t);
end