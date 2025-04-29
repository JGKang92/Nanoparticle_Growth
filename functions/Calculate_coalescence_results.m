function [sigmat_value, coal_propensity, k_f_result, k_f_exp, AtA0_theory] = Calculate_coalescence_results(x, data, rmax)
    % x(1): r_c
    % x(2): kappa_c
    % x(3): a1s (variance of S(r), nu)
    % x(4): l (ligand length)
    % x(5): P0 value
    % x(7): gamma value
    tmax = 180;
    rc = x(1);
    kappa_c = x(2);
    a1s = x(3);
    l = x(4);
    xlen = 5000;
    lent = 5000;
    tmin = 8;
    P0 = x(5);
    P0tval = x(6);
    b1k = x(7);
    timedata = data.timedata;
    AtA0_exp = data.AtA0_exp;
    sigmat = @(t) 1.7773331025798687 + 0.04422044832970424.*t - 0.0008212384775235609.*t.*t + ...
        7.651821072767395.*(10.^-6).*t.*t.*t - 2.153266857581242.*(10.^-8).*t.*t.*t.*t;
    dsigmat = @(t) 0.04422044832970424 - 0.0016424769550471218.*t + 0.000022955463218302186.*t.*t - ...
        8.613067430324968.*(10.^-8).*t.*t.*t;
    sigmat_value = sigmat(timedata);
    % clk: sigma_c, most probable core-core contact distance
    c1k = 3.007832240727688;
    % sigma-dependent coalescence propensity
    kappa_t = @(kappa_c, t) kappa_c./(exp(b1k.*(sigmat(t) - c1k))+1);
    coal_propensity = [sigmat(timedata), kappa_t(kappa_c, timedata)];
    % Reaction sink funtion
    Sr = @(r,t) exp(-(r-(sigmat(t)+l)).^2/(2*a1s));
    Pinit = @(rc,r) (r >= rc & r <= 270);
    % Solve PDE
    pdefun = @(x,t,u,dudx) coalescence_pde(x, t, u, dudx, @Dr, sigmat, dsigmat, kappa_t, Sr, kappa_c);
    icfun = @(x) coalescence_initial_cond(x, Pinit, sigmat, rc);
    bcfun = @(xl,ul,xr,ur,t) coalescence_bd_cond(xl,ul,xr,ur,t,sigmat);
    xmesh = linspace(0,rmax,xlen);
    tspan = linspace(tmin,tmax,lent);
    sol = pdepe(0,pdefun,icfun,bcfun,xmesh,tspan);
    %
    kappa_fval = zeros(1,lent);
    for idx=1:lent
        tval = tspan(idx);
        ufun_pp = csapi(xmesh, sol(idx,:));
        ufun = @(x) 4*pi*(x+sigmat(tval)).^2.*kappa_t(kappa_c,tval).*Sr(x+sigmat(tval),sigmat(tval)).*ppval(ufun_pp,x);
        kappa_fval(idx) = integral(ufun,xmesh(1),xmesh(end));
    end
    kappa_f_pp = csapi(tspan, kappa_fval);
    kappa_f = @(t) ppval(kappa_f_pp,t);
    % Coalescence rate coefficient
    k_f_result = [tspan', kappa_f(tspan)'];
    %
    kinetic_eq = @(t,y) P0.*(t< P0tval) - kappa_f(t).*y.*y;
    yinit = AtA0_exp(1);
    [~, AtA0_theory] = ode45(kinetic_eq,tspan,yinit);
    AtA0_theory = [tspan', AtA0_theory];
    % Experimental kf(t)
    tot_radi_mean_pp = csapi(8:180, AtA0_exp);
    tot_radi_mean_f = @(t) ppval(tot_radi_mean_pp,t);
    d_tot_radi_mean_pp = fnder(tot_radi_mean_pp,1);
    d_tot_radi_mean_f = @(t) fnval(d_tot_radi_mean_pp,t);
    kf_fun = @(t) (P0.*(t< P0tval) - d_tot_radi_mean_f(t))./tot_radi_mean_f(t).^2;
    k_f_exp = [(8:180)', movmean(kf_fun(8:180),20)'];
end