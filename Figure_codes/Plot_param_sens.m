function Plot_param_sens(expdata, result, a)
    %This function plots sensitivity of parameters to relative variance
    if isfield(expdata, "varnumberdata")
        relvar = expdata.varnumberdata ./ expdata.meannumberdata.^2;
    else
        relvar = expdata.varr3data./expdata.meanr3data.^2;
    end
    for idx = 1:5
        if idx == 4
            fidx = 5;
        elseif idx == 5
            fidx = 4;
        else
            fidx = idx;
        end
    
        nexttile(a + fidx)
        T = expdata.timedata;
        plot(T, relvar, ...
            'ko', 'LineWidth', 1)
        hold on
        if isfield(expdata,"tval")
            T = expdata.tval;
        end
        plot(T, result.var_theory ./ result.mean_theory.^2, ...
            '-', LineWidth =  2, Color = "#df1012")
        plot(T, result.var_theory_half{idx} ./ result.mean_theory_half{idx}.^2, ...
            '-', LineWidth = 2, Color = "#2b5f8f")
        plot(T, result.var_theory_twice{idx} ./ result.mean_theory_twice{idx}.^2, ...
            '-', LineWidth = 2, Color = "#f4bc1f")
        hold off
        switch idx
            case 1
                txt = "\Delta\epsilon_f";
            case 2
                txt = "\Delta\epsilon_e";
            case 3
                txt = "\alpha' = 4+\alpha";
            case 4
                txt = "\kappa_a\rho_{1, \infty}";
            case 5
                txt = "\kappa_a\sigma_s/D_1";
        end
        title(txt)
        
    end

end