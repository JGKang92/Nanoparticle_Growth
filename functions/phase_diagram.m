function [lambda_arr, sigmac_arr, phase_idx] = phase_diagram(rc, len)
kappa_c = 11.6;
lambda_arr = 0.5 + 1.5*rand(len,1);
sigmac_arr = 1.5 + 4.5*rand(len,1);
trange = zeros(len, 2);
parfor idx=1:len
        trange(idx,:) = Determine_trange(rc, kappa_c, lambda_arr(idx), sigmac_arr(idx));
        disp(idx)
end
phase_range = [8, 28, 41, 107, 155];
phase_idx = zeros(len, 2);
for idx = 1:len
    for idx1 = 1:(length(phase_range)-1)
        if (trange(idx, 1) >= phase_range(idx1) && trange(idx, 1) < phase_range(idx1 + 1))
            phase_idx(idx, 1) = idx1;
        end
        if (trange(idx, 2) >= phase_range(idx1) && trange(idx, 2) < phase_range(idx1 + 1))
            phase_idx(idx, 2) = idx1;
        end
    end
end
