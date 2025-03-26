function plot_curve(axis, x, y, xlbl, ylbl)
    plot(axis, y, 'k-', 'LineWidth', 1.5); 
    xlabel(axis, xlbl);
    ylabel(axis, ylbl);
    axis.Box = 'on';
    xlim(axis, [min(x), max(x)]);
    ylim(axis, [min(y)*0.9 max(y)*1.01]);
end