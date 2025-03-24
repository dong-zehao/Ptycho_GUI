function plot_probe(app, probe)
    Nmode = size(probe,2);
    Npixel = size(probe{1},1);
    % initialize the probe matrix
    probe_mat = zeros(Npixel, Npixel*Nmode);
    for imode = 1:Nmode
        idx_mode = (imode-1)*Npixel + 1 : imode*Npixel;
        probe_mat(:, idx_mode) = probe{imode};
    end
    % plot amplitude and phase respectively
    utils.imshow(app.UIAxes_3, abs(probe_mat), 'parula');
end