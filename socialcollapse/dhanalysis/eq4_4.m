% Converted from Python to MATLAB by ChatGPT
function eq4_4
    % Parameters (initial values)
    init_r = 1;
    init_alpha = 0.001;
    init_beta = 1;
    init_mu = 2;
    init_k = 6;
    init_c = 1;
    init_rho = 1;

    max_h = 10;
    max_x = 8;

    % Create figure
    f = figure('Position', [100, 100, 900, 600]);
    ax = axes('Parent', f, 'Position', [0.38, 0.1, 0.57, 0.85]);
    xlabel(ax, 'h');
    ylabel(ax, 'x');
    title(ax, sprintf('Stream: mu=%.2f, r=%.2f, alpha=%.2f, beta=%.2f, c=%.2f, k=%.2f, rho=%.2f', ...
        init_mu, init_r, init_alpha, init_beta, init_c, init_k, init_rho));

    % Grid for streamplot
    [h_grid, x_grid] = meshgrid(linspace(0, max_h, 100), linspace(0, max_x, 100));

    % Initial plot
    [dh, dx] = compute_field(h_grid, x_grid, init_r, init_alpha, init_beta, init_mu, init_k, init_c, init_rho);
    streams = streamslice(ax, h_grid, x_grid, dh, dx, 2);
    xlim(ax, [-0.1, max_h]);
    ylim(ax, [-0.1, max_x]);

    % Slider positions and ranges
    slider_params = {
        'mu (log)',    -2, 2, log10(init_mu),   0.07, 0.85, @(v)10^v;
        'r (log)',     -2, 2, log10(init_r),    0.07, 0.73, @(v)10^v;
        'alpha (log)', -3, 0, log10(init_alpha),0.07, 0.61, @(v)10^v;
        'beta (log)',   0, 3, log10(init_beta), 0.07, 0.49, @(v)10^v;
        'c (log)',     -1, 1, log10(init_c),    0.07, 0.37, @(v)10^v;
        'k (log)',      0, 1, log10(init_k),    0.07, 0.25, @(v)10^v;
        'rho (log)',   -1, 1, log10(init_rho),  0.07, 0.13, @(v)10^v;
    };

    sliders = gobjects(1,7);
    for i = 1:7
        sliders(i) = uicontrol('Style', 'slider', ...
            'Min', slider_params{i,2}, 'Max', slider_params{i,3}, ...
            'Value', slider_params{i,4}, ...
            'Units', 'normalized', ...
            'Position', [slider_params{i,5}, slider_params{i,6}, 0.18, 0.07], ...
            'Callback', @(src,~) update_plot());
        uicontrol('Style', 'text', ...
            'Units', 'normalized', ...
            'Position', [slider_params{i,5}, slider_params{i,6}+0.06, 0.18, 0.03], ...
            'String', slider_params{i,1});
    end

    function update_plot()
        mu    = slider_params{1,7}(get(sliders(1), 'Value'));
        r     = slider_params{2,7}(get(sliders(2), 'Value'));
        alpha = slider_params{3,7}(get(sliders(3), 'Value'));
        beta  = slider_params{4,7}(get(sliders(4), 'Value'));
        c     = slider_params{5,7}(get(sliders(5), 'Value'));
        k     = slider_params{6,7}(get(sliders(6), 'Value'));
        rho   = slider_params{7,7}(get(sliders(7), 'Value'));

        [dh, dx] = compute_field(h_grid, x_grid, r, alpha, beta, mu, k, c, rho);
        cla(ax);
        streamslice(ax, h_grid, x_grid, dh, dx, 2);
        xlim(ax, [-0.1, max_h]);
        ylim(ax, [-0.1, max_x]);
        xlabel(ax, 'h');
        ylabel(ax, 'x');
        title(ax, sprintf('Stream: mu=%.2f, r=%.2f, alpha=%.2f, beta=%.2f, c=%.2f, k=%.2f, rho=%.2f', ...
            mu, r, alpha, beta, c, k, rho));
    end
end

function [dh, dx] = compute_field(h, x, r, alpha, beta, mu, k, c, rho)
    gamma = @(x) beta - (beta - alpha) .* x ./ k;
    dh = h .* (r - gamma(x) .* h);
    dx = mu .* x .* (1 - x ./ k) - c .* x .* h ./ (rho + x);
end