function [result, fitdata, expdata] = post_processing(result, fitdata, expdata, options, ...
    smoothing_indices, overwrite_indices, overwrite_values, situ, TF, t0idx)
    if ~isvector(fitdata.xfit)
        error("xfit should be vector")
    end
    % smoothing
    for ii = 1:numel(smoothing_indices)
        indices = smoothing_indices{ii};
        fitdata.xfit(indices) = smoothdata(fitdata.xfit(indices),'gaussian', 5);
    end
    
    % Overwriting
    for ii = 1:numel(overwrite_indices)
        fitdata.xfit(overwrite_indices(ii)) = overwrite_values(ii);
    end
    
    % Calculate theoretical results
    order_max = 10;
    lent = length(expdata.timedata);
    
    switch situ
        case "in"
            % Interpolate parameters
            nIndices = cellfun(@numel, options.indices);
            lent_ = unique(nIndices(nIndices~=1));
            if ~isscalar(lent_)
                error("lent should be scalar")
            end
            if numel(expdata.timedata) ~= lent_
                if isfield(options,"bs")
                    t_old = expdata.timedata(1:bs:end);
                else
                    t_old = expdata.timedata(1:2:end);
                end
                [xfit_new, ind_new, lb_new, ub_new] = fill_parameters(fitdata.xfit, options.indices,t_old, expdata.timedata, fitdata.lb, fitdata.ub, "linear","linear");
                options.indices = ind_new;
                fitdata.other_opt = options;
                fitdata.xfit = xfit_new;
                fitdata.lb = lb_new;
                fitdata.ub = ub_new;
            end
            
            [result.mean_theory, result.var_theory, result.totmdata] = ...
                Calculate_Mean_Var_In_situ(fitdata.xfit, expdata, options);
            disp("Calculated Mean, Variance, and total monomer concentration")

            [result.Jntilde, result.psi_n_val, result.rhoct, result.ratio_comparison] = ...
                Calculate_Jn_tilde_In_situ(fitdata.xfit, expdata, result, options);
            disp("Calculated Jn tilde")
            
            result.moments = zeros(lent, order_max);
            result.moments(:,1) = result.mean_theory;
            result.moments(:,2) = result.var_theory + result.mean_theory.^2;
            expdata.moment_raw_value = zeros(lent, order_max);
            expdata.moment_raw_value(:,1) = expdata.meannumberdata;
            expdata.moment_raw_value(:,2) = expdata.varnumberdata + expdata.meannumberdata.^2;
            for order = 3:order_max
                [expdata.moment_raw_value(:,order), result.moments(:,order)] = ...
                    Calculate_higher_moment_In_situ(fitdata.xfit, expdata, order, options);
            end
            disp("Calculated higher moments")
            
            % Using higher moments, predict size distribution
            % Note: This part works better in Mathematica.
            result.pnt_guess = Construct_pnt(expdata);
            disp("Calculated pnt using higher moments")

            % Obtain Jntilde using above pnt
            [result.Jntilde_pnt_guess, ~, ~, ~] = ...
                Calculate_Jn_tilde_In_situ(fitdata.xfit, expdata, result, options);
            disp("Calculated Jn tilde")

            [result.mun_s, result.monomer_chemical_potential, result.free_energy] = ...
                Calculate_chemical_potential(fitdata.xfit, expdata, options);
            disp("Calculated chemical potential")
            
            result.pw_growth_rate = Calculate_pnjn(fitdata.xfit, ...
                expdata, result, options);
            disp("Calculated pnjn")

            [result.exp_mun, result.muns_mn, result.muns_log] = ...
                Extract_chemical_potential_In_situ(fitdata.xfit, expdata, options);
            disp("Extracted chemical potential")

            result = Calculate_param_sens("in",expdata, fitdata, result);
            disp("Calculated Parameter sensitivity")

            if TF(1)
                % Cuboctahedron case, shape_index = 1
                indices = options.indices;
                xfit = fitdata.xfit;
                fitdata.xfit_modified = xfit;
                fitdata.xfit_modified(indices{8}) = ones(size(indices{8}));
                [result.mean_theory_cOh, result.var_theory_cOh, ~] = ...
                    Calculate_Mean_Var_In_situ(fitdata.xfit_modified, expdata, options);
                disp("Calculated mean and variance assuming shape_index = 1")

                % Truncated octahedron case, shape_index = 4/3
                fitdata.xfit_modified = xfit;
                fitdata.xfit_modified(indices{8}) = ones(size(indices{8}))*4/3;    
                [result.mean_theory_tOh, result.var_theory_tOh, ~] = ...
                     Calculate_Mean_Var_In_situ(fitdata.xfit_modified, expdata, options);
                disp("Calculated mean and variance assuming shape_index = 4/3")

                % Truncated octahedron case, shape_index = 1.8
                fitdata.xfit_modified = xfit;
                fitdata.xfit_modified(indices{8}) = ones(size(indices{8}))*1.8;    
                [result.mean_theory_18, result.var_theory_18, ~] = ...
                    Calculate_Mean_Var_In_situ(fitdata.xfit_modified, expdata, options);
                disp("Calculated mean and variance assuming shape_index = 1.8")

                % Octahedron case, shape_index = 2
                fitdata.xfit_modified = xfit;
                fitdata.xfit_modified(indices{8}) = ones(size(indices{8}))*2;
                [result.mean_theory_Oh, result.var_theory_Oh, ~] = ...
                     Calculate_Mean_Var_In_situ(fitdata.xfit_modified, expdata, options);
                
                disp("Calculated mean and variance assuming shape_index = 2")
            end
        case "ex2"
            sigma_s = fitdata.xfit(fitdata.other_opt.indices{6});
            expdata = r2n_optim(expdata, sigma_s, fitdata.other_opt.IP_method, options.smooth_window);
            disp("Changed rdata to ndata")
            
            [result.mean_theory, result.var_theory, result.totmdata] = ...
                Calculate_Mean_Var_Ex_situ2(fitdata.xfit, expdata, options);
            disp("Calculated Mean, Variance, and total monomer concentration")

            [result.mun_s, result.monomer_chemical_potential, result.free_energy] = ...
                Calculate_chemical_potential(fitdata.xfit, expdata, options);
            disp("Calculated chemical potential")
            
            [result.exp_mun, result.muns_mn, result.muns_log] = ...
                Extract_chemical_potential_Ex_situ2(fitdata.xfit, expdata, options);
            disp("Extracted chemical potential")

            [result.Jntilde, result.psi_n_val, result.rhoct] = ...
                Calculate_Jn_tilde_Ex_situ2(fitdata.xfit, expdata, result, options);
            disp("Calculated Jn tilde")

            expdata.moment_raw_value = zeros(lent, order_max);
            expdata.moment_raw_value(:,1) = expdata.meanr3data./sigma_s^3;
            expdata.moment_raw_value(:,2) = (expdata.varr3data + expdata.meanr3data.^2)./sigma_s^6;
            for idx = 3:order_max
                mmt = Calculate_higher_moment_Ex_situ2(fitdata.xfit, expdata, idx, options);
                expdata.moment_raw_value(:,idx) = mmt;
            end
            disp("Calculated higher moments")

            % Predict size distribution using higher moments
            result.pnt_guess = Construct_pnt(expdata);
            disp("Calculated pnt using higher moments")

            result = Calculate_param_sens("ex2",expdata, fitdata, result);
            disp("Calculated Parameter sensitivity")
        otherwise
            error("wrong input")
    end
    if ~isnan(t0idx)
        [expdata.relvar, result.eq_m12_relvar, result.eta2_inf] = Calculate_approximate_relvar(expdata, result, t0idx);
    end
    disp("Finished post-processing")
end