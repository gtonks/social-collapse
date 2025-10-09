% Analytic Version: Uses fsolve for equilibria
function Social_Collapse
    % Parameters
    init_r = 2.304;
    init_beta = 30;
    init_alpha = 2.4;
    init_c = 0.600;
    init_rho = 0.02;
    init_max_H = 1.5;
    init_max_x = 1.180;

    % Create figures
    f1 = figure('Position', [100, 100, 600, 600]);
    ax = axes('Parent', f1, 'Position', [0.1, 0.1, 0.85, 0.85]);
    xlabel(ax, 'Resources, y');
    ylabel(ax, 'Population size, H');

    f2 = figure('Position', [800, 100, 600, 600]);

    % Grid
    [x_grid, H_grid] = meshgrid(linspace(0, init_max_x, 100), linspace(0, init_max_H, 100));

    % Initial plot
    [dH, dx] = compute_field(H_grid, x_grid, init_r, init_alpha, init_beta, init_c, init_rho);
    streamslice(ax, x_grid, H_grid, dx, dH, 2, 'Color', 'r');
    hold(ax, 'on');
    [eq_points, is_stable, eq_types] = compute_equilibria(init_r, init_alpha, init_beta, init_c, init_rho, init_max_x, init_max_H);
    valid_idx = ~isnan(eq_points(:,1)) & ~isinf(eq_points(:,1)) & ...
                eq_points(:,1) >= 0 & eq_points(:,1) <= init_max_x & ...
                ~isnan(eq_points(:,2)) & ~isinf(eq_points(:,2)) & ...
                eq_points(:,2) >= 0 & eq_points(:,2) <= init_max_H;
    plot(ax, eq_points(valid_idx & is_stable, 1), eq_points(valid_idx & is_stable, 2), ...
         'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot(ax, eq_points(valid_idx & ~is_stable, 1), eq_points(valid_idx & ~is_stable, 2), ...
         'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
    hold(ax, 'off');
    xlim(ax, [-0.1, init_max_x]);
    ylim(ax, [-0.1, init_max_H]);

    % Equilibria list
    eq_text = uicontrol('Parent', f2, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.05, 0.05, 0.45, 0.90], ...
        'String', format_equilibria(eq_points, is_stable, eq_types, valid_idx), ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'white');

    % Sliders
    slider_params = {
        'r',     0, 10,   init_r,     0.55, 0.80, @(v)v, 'r';
        'alpha', 0, init_beta, init_alpha, 0.55, 0.68, @(v)v, 'a';
        'beta',  0, 60,   init_beta,  0.55, 0.56, @(v)v, 'b';
        'c',     0, 2,    init_c,     0.55, 0.44, @(v)v, 'C';
        'rho',   0, 0.1,  init_rho,   0.55, 0.32, @(v)v, 'rho';
        'max_x', 0.5, 1.5, init_max_x, 0.55, 0.20, @(v)v, 'max_x';
        'max_H', 0.1, 1.5,  init_max_H, 0.55, 0.08, @(v)v, 'max_H';
    };

    sliders = gobjects(1,7);
    value_labels = gobjects(1,7);
    range_labels = gobjects(1,7);
    name_labels = gobjects(1,7);
    for i = 1:7
        sliders(i) = uicontrol('Parent', f2, 'Style', 'slider', ...
            'Min', slider_params{i,2}, 'Max', slider_params{i,3}, ...
            'Value', slider_params{i,4}, ...
            'Units', 'normalized', ...
            'Position', [slider_params{i,5}, slider_params{i,6}, 0.30, 0.05], ...
            'Callback', @(src,~) update_equilibria());
        name_labels(i) = uicontrol('Parent', f2, 'Style', 'text', ...
            'Units', 'normalized', ...
            'Position', [slider_params{i,5}-0.15, slider_params{i,6}+0.06, 0.10, 0.03], ...
            'String', slider_params{i,8}, ...
            'HorizontalAlignment', 'right');
        value_labels(i) = uicontrol('Parent', f2, 'Style', 'text', ...
            'Units', 'normalized', ...
            'Position', [slider_params{i,5}+0.31, slider_params{i,6}, 0.09, 0.05], ...
            'String', sprintf('%.3f', slider_params{i,4}), ...
            'HorizontalAlignment', 'left');
        range_labels(i) = uicontrol('Parent', f2, 'Style', 'text', ...
            'Units', 'normalized', ...
            'Position', [slider_params{i,5}-0.04, slider_params{i,6}+0.06, 0.10, 0.03], ...
            'String', sprintf('[%.2f, %.2f]', slider_params{i,2}, slider_params{i,3}), ...
            'HorizontalAlignment', 'left');
    end

    function update_equilibria()
        r     = slider_params{1,7}(get(sliders(1), 'Value'));
        alpha = slider_params{2,7}(get(sliders(2), 'Value'));
        beta  = slider_params{3,7}(get(sliders(3), 'Value'));
        c     = slider_params{4,7}(get(sliders(4), 'Value'));
        rho   = slider_params{5,7}(get(sliders(5), 'Value'));
        max_x = slider_params{6,7}(get(sliders(6), 'Value'));
        max_H = slider_params{7,7}(get(sliders(7), 'Value'));

        set(sliders(2), 'Max', beta);
        if alpha > beta
            set(sliders(2), 'Value', beta);
            alpha = beta;
        end

        for i = 1:7
            set(value_labels(i), 'String', sprintf('%.3f', slider_params{i,7}(get(sliders(i), 'Value'))));
        end
        set(range_labels(2), 'String', sprintf('[%.2f, %.2f]', 0, beta));

        [x_grid, H_grid] = meshgrid(linspace(0, max_x, 100), linspace(0, max_H, 100));

        cla(ax);
        [dH, dx] = compute_field(H_grid, x_grid, r, alpha, beta, c, rho);
        streamslice(ax, x_grid, H_grid, dx, dH, 2, 'Color', 'r');
        hold(ax, 'on');
        [eq_points, is_stable, eq_types] = compute_equilibria(r, alpha, beta, c, rho, max_x, max_H);
        valid_idx = ~isnan(eq_points(:,1)) & ~isinf(eq_points(:,1)) & ...
                    eq_points(:,1) >= 0 & eq_points(:,1) <= max_x & ...
                    ~isnan(eq_points(:,2)) & ~isinf(eq_points(:,2)) & ...
                    eq_points(:,2) >= 0 & eq_points(:,2) <= max_H;
        plot(ax, eq_points(valid_idx & is_stable, 1), eq_points(valid_idx & is_stable, 2), ...
             'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
        plot(ax, eq_points(valid_idx & ~is_stable, 1), eq_points(valid_idx & ~is_stable, 2), ...
             'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
        hold(ax, 'off');
        xlim(ax, [-0.1, max_x]);
        ylim(ax, [-0.1, max_H]);
        xlabel(ax, 'Resources, x');
        ylabel(ax, 'Population size, H');

        set(eq_text, 'String', format_equilibria(eq_points, is_stable, eq_types, valid_idx));
    end
end

function [dH, dx] = compute_field(H, x, r, alpha, beta, c, rho)
    gamma = @(x) beta - (beta - alpha) .* x;
    dH = H .* (r - gamma(x) .* H);
    dx = x .* (1 - x) - c .* x .* H ./ (rho + x);
end

function [eq_points, is_stable, eq_types] = compute_equilibria(r, alpha, beta, c, rho, max_x, max_H)
    fprintf('Parameters: r=%.3f, alpha=%.3f, beta=%.3f, c=%.3f, rho=%.3f, max_x=%.3f, max_H=%.3f\n', ...
            r, alpha, beta, c, rho, max_x, max_H);
    
    eq_points = NaN(6, 2);
    is_stable = false(6, 1);
    eq_types = cell(6, 1);
    
    % Trivial equilibria
    eq_points(1, :) = [0, 0];
    eq_points(2, :) = [1, 0];
    eq_points(3, :) = [0, r / beta];
    
    % Stability for (0,0)
    J = [r, 0; 0, 1];
    eigvals = [r; 1];
    is_stable(1) = all(real(eigvals) < -1e-10);
    eq_types{1} = classify_equilibrium(eigvals, is_stable(1));
    fprintf('Equilibrium 1 (x=0, H=0): Jacobian = [%.12f, %.12f; %.12f, %.12f], Eigenvalues = [%.12f, %.12f], Stable=%d, Type=%s\n', ...
            J(1,1), J(1,2), J(2,1), J(2,2), eigvals(1), eigvals(2), is_stable(1), eq_types{1});
    
    % Stability for (1,0)
    J = [r, 0; -c / (rho + 1), -1]; % Corrected Jacobian
    eigvals = [r; -1];
    is_stable(2) = all(real(eigvals) < -1e-10);
    eq_types{2} = classify_equilibrium(eigvals, is_stable(2));
    fprintf('Equilibrium 2 (x=1, H=0): Jacobian = [%.12f, %.12f; %.12f, %.12f], Eigenvalues = [%.12f, %.12f], Stable=%d, Type=%s\n', ...
            J(1,1), J(1,2), J(2,1), J(2,2), eigvals(1), eigvals(2), is_stable(2), eq_types{2});
    
    % Stability for (0, r/beta)
    H_eq = r / beta;
    J = [-r, (beta - alpha) * (r / beta)^2; 0, 1 - (c * r) / (beta * rho)];
    eigvals = [-r; 1 - (c * r) / (beta * rho)];
    is_stable(3) = all(real(eigvals) < -1e-10);
    eq_types{3} = classify_equilibrium(eigvals, is_stable(3));
    fprintf('Equilibrium 3 (x=0, H=%.12f): Jacobian = [%.12f, %.12f; %.12f, %.12f], Eigenvalues = [%.12f, %.12f], Stable=%d, Type=%s\n', ...
            H_eq, J(1,1), J(1,2), J(2,1), J(2,2), eigvals(1), eigvals(2), is_stable(3), eq_types{3});
    
    % Non-trivial equilibria
    ode_fun = @(x) x * (1 - x) - c * x * (r / (beta - (beta - alpha) * x)) / (rho + x);
    x_guesses = [linspace(0.01, 0.99, 500), linspace(0.05, 0.15, 200), linspace(0.6, 0.7, 200)];
    
    fprintf('Debug: Initializing non-trivial equilibria search\n');
    x_sols = [];
    H_sols = [];
    for x0 = x_guesses
        try
            x_sol = fsolve(ode_fun, x0, optimoptions('fsolve', 'Display', 'none', ...
                'TolFun', 1e-14, 'TolX', 1e-14, 'MaxIterations', 5000));
            if ~isnan(x_sol) && ~isinf(x_sol) && x_sol > 1e-8 && x_sol < 1 - 1e-8
                gamma_x = beta - (beta - alpha) * x_sol;
                if gamma_x > 1e-8
                    H_sol = r / gamma_x;
                    if ~isnan(H_sol) && ~isinf(H_sol) && H_sol > 1e-8 && H_sol <= max_H
                        residual = abs(ode_fun(x_sol));
                        if residual < 1e-10
                            if (abs(x_sol) > 1e-6 || abs(H_sol - r/beta) > 1e-6) && ...
                               (isempty(x_sols) || min(abs(x_sols - x_sol)) > 1e-6)
                                x_sols = [x_sols; x_sol];
                                H_sols = [H_sols; H_sol];
                                fprintf('Root candidate: x=%.12f, H=%.12f, gamma_x=%.12f, residual=%.2e\n', ...
                                        x_sol, H_sol, gamma_x, residual);
                            end
                        end
                    end
                end
            end
        catch
            continue;
        end
    end
    [x_sols, idx] = sort(x_sols);
    H_sols = H_sols(idx);
    fprintf('Found %d non trivial solutions:\n', length(x_sols));
    for i = 1:length(x_sols)
        fprintf('Solution %d: x=%.12f, H=%.12f\n', i, x_sols(i), H_sols(i));
    end
    
    % Store non-trivial
    for i = 1:min(length(x_sols), 3)
        eq_points(3+i, :) = [x_sols(i), H_sols(i)];
        x_eq = x_sols(i);
        H_eq = H_sols(i);
        gamma_x = beta - (beta - alpha) * x_eq;
        J = [r - 2 * gamma_x * H_eq, H_eq^2 * (beta - alpha); ...
             -c * x_eq / (rho + x_eq), 1 - 2 * x_eq - c * H_eq * rho / (rho + x_eq)^2];
        fprintf('Jacobian at (x=%.12f, H=%.12f):\n', x_eq, H_eq);
        fprintf('[%.12f, %.12f]\n[%.12f, %.12f]\n', J(1,1), J(1,2), J(2,1), J(2,2));
        eigvals = eig(J);
        is_stable(3+i) = all(real(eigvals) < -1e-10);
        eq_types{3+i} = classify_equilibrium(eigvals, is_stable(3+i));
        fprintf('Equilibrium %d (x=%.12f, H=%.12f): Eigenvalues = [%.12f, %.12f], Stable=%d, Type=%s\n', ...
                3+i, x_eq, H_eq, eigvals(1), eigvals(2), is_stable(3+i), eq_types{3+i});
    end
end

function type = classify_equilibrium(eigvals, is_stable)
    if any(abs(imag(eigvals)) > 1e-10)
        if is_stable
            type = 'Stable Spiral';
        else
            type = 'Unstable Spiral';
        end
    else
        if eigvals(1) * eigvals(2) < -1e-10
            type = 'Saddle Point';
        else
            if is_stable
                type = 'Stable Node';
            else
                type = 'Unstable Node';
            end
        end
    end
end

function str = format_equilibria(eq_points, is_stable, eq_types, valid_idx)
    str = {'Equilibria'};
    count = 1;
    for i = 1:size(eq_points, 1)
        if valid_idx(i)
            x = eq_points(i, 1);
            H = eq_points(i, 2);
            stability = 'Unstable';
            if is_stable(i)
                stability = 'Stable';
            end
            str{end+1} = sprintf('%d. (%.6f, %.6f) - %s, %s', count, x, H, stability, eq_types{i});
            count = count + 1;
        end
    end
    if count == 1
        str{end+1} = '';
    end
end