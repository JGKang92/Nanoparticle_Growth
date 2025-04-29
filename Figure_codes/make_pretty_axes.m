function [ax, xl, yl] = make_pretty_axes(ax)
    if isequal(ax.XColor, ones(1, 3)*0.15) % default setting of XColor
        set(ax, XColor = zeros(1, 3)) % Change XColor to black
    end
    if isequal(ax.YColor, ones(1, 3)*0.15) % default setting of YColor
        set(ax, YColor = zeros(1, 3)) % Change YColor to black
    end
    set(ax, Color = 'none', ...
            LineWidth = 1, Box = "off", ...
            FontSize = 16, ...
            TickLength = [0.02, 0.01])
     
    if strcmp(ax.XAxisLocation,"bottom")
        xl = xlim(ax);
        if isinf(xl(end))
            TF = arrayfun(@(x)isprop(x,"XData"),ax.Children);
            x_end = max(arrayfun(@(x)max(x.XData),ax.Children(TF)));
        else
            x_end = xl(end);
        end
        xline(ax, x_end, LineWidth = ax.LineWidth, Color = ax.YColor, Alpha = 1, DisplayName = "")
    end
    if strcmp(ax.YAxisLocation, "left")
        yl = ylim(ax);
        if isinf(yl(end))
            TF = arrayfun(@(y)isprop(y,"YData"),ax.Children);
            y_end = max(arrayfun(@(y)max(y.YData), ax.Children(TF)));
        else
            y_end = yl(end);
        end
        yline(ax, y_end, LineWidth = ax.LineWidth, Color = ax.XColor, Alpha = 1, DisplayName = "")
    end
end