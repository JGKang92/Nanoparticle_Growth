function Plot_higher_moments(ax, expdata, result, order_x, skipPoint)
    T = expdata.timedata;
    for idx = 3:order_x
        p = semilogy(ax, T(1:skipPoint:end), expdata.moment_raw_value(1:skipPoint:end,idx), 'o', ...
            'MarkerSize', 10);
        hold on
        semilogy(ax, T, result.moments(:,idx), '-', ...
            'LineWidth', 2, ...
            'Color',p.Color)
    end
    hold off
end