function Plot_Moment(ax, ed, fd, situ, name, color)
    t = ed.timedata;
    switch situ
        case "in"
            if ~isvector(fd.xfit)
                error("fd.xfit should be vector")
            end
            [mean_theory, var_theory] = Calculate_Mean_Var_In_situ(fd.xfit, ed, fd.other_opt);
        case "ex2"
    end
    factor = fd.other_opt.factor;
    switch name
        case {"mean","first_moment"}
            y = ed.meannumberdata;
            y1 = ed.smooth_meannumberdata;
            y2 = mean_theory;
        case {"variance","var"}
            y = ed.varnumberdata;
            y1 = factor*ed.smooth_varnumberdata;
            y2 = var_theory;
        case {"relative_variance", "rv","relvar"}
            y = ed.varnumberdata./ed.meannumberdata.^2;
            y1 = factor*ed.smooth_varnumberdata./ed.smooth_meannumberdata.^2;
            y2 = var_theory./mean_theory.^2;
        case {"Coefficient_of_Variation","CV"}
            y = sqrt(ed.varnumberdata)./ed.meannumberdata;
            y1 = sqrt(factor*ed.smooth_varnumberdata)./ed.smooth_meannumberdata;
            y2 = sqrt(var_theory)./mean_theory;
    end
    plot(ax, t, y,'s', Color = 0.8.*ones(1,3))
    ax.NextPlot = "add";
    plot(ax, t, y1,'ko')
    plot(ax, t, y2, Color = color, LineWidth = 2)
    title(name)
end