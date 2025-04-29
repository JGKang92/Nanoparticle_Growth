function result_out = Calculate_param_sens(situ, expdata, fitdata, result)
    %This function plots sensitivity of parameters to relative variance
    result_out = result;
    mean_theory_half = cell(1, 5);
    var_theory_half = cell(1, 5);
    mean_theory_twice = cell(1, 5);
    var_theory_twice = cell(1, 5);
    for idx = 1:5
        xfit_halved = fitdata.xfit;
        xfit_halved(idx) = xfit_halved(idx)./2;
        switch situ
            case "in"
                [mean_theory_half{idx}, var_theory_half{idx}, ~] = ...
                    Calculate_Mean_Var_In_situ(xfit_halved, expdata, fitdata.other_opt);
            case "ex2"
                [mean_theory_half{idx}, var_theory_half{idx}, ~] = ...
                    Calculate_Mean_Var_Ex_situ2(xfit_halved, expdata, fitdata.other_opt);
        end
        
        xfit_twiced = fitdata.xfit;
        xfit_twiced(idx) = xfit_twiced(idx)*2;
        switch situ
            case "in"
                [mean_theory_twice{idx}, var_theory_twice{idx}, ~] = ...
                    Calculate_Mean_Var_In_situ(xfit_twiced, expdata, fitdata.other_opt);
            case "ex2"
                [mean_theory_twice{idx}, var_theory_twice{idx}, ~] = ...
                    Calculate_Mean_Var_Ex_situ2(xfit_twiced, expdata, fitdata.other_opt);
        end
    end
    result_out.mean_theory_half = mean_theory_half;
    result_out.var_theory_half = var_theory_half;
    result_out.mean_theory_twice = mean_theory_twice;
    result_out.var_theory_twice = var_theory_twice;
end