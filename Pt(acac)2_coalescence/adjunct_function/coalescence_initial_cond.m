function u0 = coalescence_initial_cond(x, Pinit, sigmat, sigmac)
    % Pinit0: Initial distribution, normalized by
    % integral(Sn*pi*r^(d-1)*Pinit(r), 0, inf) = 1
    % sigmat: sigma value in function handle
    u0 = Pinit(sigmac, x + sigmat(0));
end