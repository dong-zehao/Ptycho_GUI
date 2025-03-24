function plot_error(app, fourier_error)
    x_iters=[1:size(fourier_error,1)];
    fourier_error_out = mean(gather(fourier_error),2,'omitnan');
    y_ferr = fourier_error_out(~isnan(fourier_error_out));
    utils.plot_curve(app.UIAxes, x_iters, y_ferr, 'Iteration', 'Fourier error')
end