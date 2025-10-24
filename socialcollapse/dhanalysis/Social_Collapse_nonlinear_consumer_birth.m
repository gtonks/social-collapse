% Analytic Version: Uses fsolve for equilibria
function Social_Collapse_B
    % Parameters
    init_B = 3.630;
    init_d = 0.266;
    init_beta = 12.444;
    init_alpha = 1.080;
    init_c = 0.600;
    init_rho = 0.033;
    init_max_H = 0.500;
    init_max_x = 1.0500;

    % Create figures
    f1 = figure('Position', [100, 100, 600, 600]);
    ax = axes('Parent', f1, 'Position', [0.1, 0.1, 0.85, 0.85]);
    xlabel(ax, 'Resources, y', 'FontSize',16);
    ylabel(ax, 'Consumer population, H','FontSize',16);
    xticks([0 1]);
    yticks([0 0.5])

    f2 = figure('Position', [800, 100, 600, 600]);

    % Grid
    [x_grid, H_grid] = meshgrid(linspace(0, init_max_x, 100), linspace(0, init_max_H, 100));

    % Initial plot
    [dH, dx] = compute_field(H_grid, x_grid, init_B, init_d, init_alpha, init_beta, init_c, init_rho);
    streamslice(ax, x_grid, H_grid, dx, dH, 2, 'Color', 'r');
    hold(ax, 'on');
    [eq_points, is_stable, eq_types] = compute_equilibria(init_B, init_d, init_alpha, init_beta, init_c, init_rho, init_max_x, init_max_H);
    valid_idx = ~isnan(eq_points(:,1)) & ~isinf(eq_points(:,1)) & ...
                eq_points(:,1) >= 0 & eq_points(:,1) <= init_max_x & ...
                ~isnan(eq_points(:,2)) & ~isinf(eq_points(:,2)) & ...
                eq_points(:,2) >= 0 & eq_points(:,2) <= init_max_H;
    plot(ax, eq_points(valid_idx & is_stable, 1), eq_points(valid_idx & is_stable, 2), ...
         'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot(ax, eq_points(valid_idx & ~is_stable, 1), eq_points(valid_idx & ~is_stable, 2), ...
         'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
    hold(ax, 'off');
    xlim(ax, [-0.05, init_max_x]);
    ylim(ax, [-0.01, init_max_H]);
    

    % Equilibria list
    eq_text = uicontrol('Parent', f2, 'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.05, 0.05, 0.45, 0.90], ...
        'String', format_equilibria(eq_points, is_stable, eq_types, valid_idx), ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'white');

    % Sliders
    slider_params = {
        'B',     0, 5,    init_B,     0.55, 0.80, @(v)v, 'B';
        'd',     0, 2,    init_d,     0.55, 0.68, @(v)v, 'd';
        'alpha', 0, init_beta, init_alpha, 0.55, 0.56, @(v)v, 'a';
        'beta',  0, 60,   init_beta,  0.55, 0.44, @(v)v, 'b';
        'c',     0, 2,    init_c,     0.55, 0.32, @(v)v, 'C';
        'rho',   0, 0.1,  init_rho,   0.55, 0.20, @(v)v, 'rho';
        'max_x', 0.5, 1.5, init_max_x, 0.55, 0.14, @(v)v, 'max_x';
        'max_H', 0.1, 1.5, init_max_H, 0.55, 0.08, @(v)v, 'max_H';
    };

    sliders = gobjects(1,8);
    value_labels = gobjects(1,8);
    range_labels = gobjects(1,8);
    name_labels = gobjects(1,8);
    for i = 1:8
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
        B     = slider_params{1,7}(get(sliders(1), 'Value'));
        d     = slider_params{2,7}(get(sliders(2), 'Value'));
        alpha = slider_params{3,7}(get(sliders(3), 'Value'));
        beta  = slider_params{4,7}(get(sliders(4), 'Value'));
        c     = slider_params{5,7}(get(sliders(5), 'Value'));
        rho   = slider_params{6,7}(get(sliders(6), 'Value'));
        max_x = slider_params{7,7}(get(sliders(7), 'Value'));
        max_H = slider_params{8,7}(get(sliders(8), 'Value'));

        set(sliders(3), 'Max', beta);
        if alpha > beta
            set(sliders(3), 'Value', beta);
            alpha = beta;
        end

        for i = 1:8
            set(value_labels(i), 'String', sprintf('%.3f', slider_params{i,7}(get(sliders(i), 'Value'))));
        end
        set(range_labels(3), 'String', sprintf('[%.2f, %.2f]', 0, beta));

        [x_grid, H_grid] = meshgrid(linspace(0, max_x, 100), linspace(0, max_H, 100));

        cla(ax);
        [dH, dx] = compute_field(H_grid, x_grid, B, d, alpha, beta, c, rho);
        streamslice(ax, x_grid, H_grid, dx, dH, 2, 'Color', 'r');
        hold(ax, 'on');
        [eq_points, is_stable, eq_types] = compute_equilibria(B, d, alpha, beta, c, rho, max_x, max_H);
        valid_idx = ~isnan(eq_points(:,1)) & ~isinf(eq_points(:,1)) & ...
                    eq_points(:,1) >= 0 & eq_points(:,1) <= max_x & ...
                    ~isnan(eq_points(:,2)) & ~isinf(eq_points(:,2)) & ...
                    eq_points(:,2) >= 0 & eq_points(:,2) <= max_H;
        plot(ax, eq_points(valid_idx & is_stable, 1), eq_points(valid_idx & is_stable, 2), ...
             'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
        plot(ax, eq_points(valid_idx & ~is_stable, 1), eq_points(valid_idx & ~is_stable, 2), ...
             'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
        hold(ax, 'off');
        xlim(ax, [-0.01, max_x]);
        ylim(ax, [-0.01, max_H]);
        xlabel(ax, 'Resources, y');
        ylabel(ax, 'Population size, H');

        set(eq_text, 'String', format_equilibria(eq_points, is_stable, eq_types, valid_idx));
    end
end

function [dH, dx] = compute_field(H, x, B, d, alpha, beta, c, rho)
    gamma = @(x) beta - (beta - alpha) .* x;
    dH = (B .* x .* H) ./ (rho + x) - d .* H - gamma(x) .* H.^2;
    dx = x .* (1 - x) - c .* x .* H ./ (rho + x);
end

function [eq_points, is_stable, eq_types] = compute_equilibria(B, d, alpha, beta, c, rho, max_x, max_H)
    fprintf('Parameters: B=%.3f, d=%.3f, alpha=%.3f, beta=%.3f, c=%.3f, rho=%.3f, max_x=%.3f, max_H=%.3f\n', ...
            B, d, alpha, beta, c, rho, max_x, max_H);
    
    eq_points = NaN(6, 2);
    is_stable = false(6, 1);
    eq_types = cell(6, 1);
    
    % Trivial equilibria
    eq_points(1, :) = [0, 0];
    eq_points(2, :) = [1, 0];
    
    % Stability for (0,0)
    J = [-d, 0; 0, 1];
    eigvals = [-d; 1];
    is_stable(1) = all(real(eigvals) < -1e-8);
    eq_types{1} = classify_equilibrium(eigvals, is_stable(1));
    fprintf('Equilibrium 1 (x=0, H=0): Jacobian = [%.12f, %.12f; %.12f, %.12f], Eigenvalues = [%.12f, %.12f], Stable=%d, Type=%s\n', ...
            J(1,1), J(1,2), J(2,1), J(2,2), eigvals(1), eigvals(2), is_stable(1), eq_types{1});
    
    % Stability for (1,0)
    J = [B/(rho+1) - d, 0; -c / (rho + 1), -1];
    eigvals = [B/(rho+1) - d; -1];
    is_stable(2) = all(real(eigvals) < -1e-8);
    eq_types{2} = classify_equilibrium(eigvals, is_stable(2));
    fprintf('Equilibrium 2 (x=1, H=0): Jacobian = [%.12f, %.12f; %.12f, %.12f], Eigenvalues = [%.12f, %.12f], Stable=%d, Type=%s\n', ...
            J(1,1), J(1,2), J(2,1), J(2,2), eigvals(1), eigvals(2), is_stable(2), eq_types{2});
    
    % Non-trivial equilibria
    ode_fun = @(x) x * (1 - x) - c * x / (rho + x) * ((B * x / (rho + x) - d) / (beta - (beta - alpha) * x));
    x_guesses = [linspace(0.01, 0.99, 500), linspace(0.05, 0.15, 200), linspace(0.6, 0.7, 200), linspace(0.75, 0.85, 200)];
    
    fprintf('Debug: Initializing non-trivial equilibria search\n');
    x_sols = [];
    H_sols = [];
    for x0 = x_guesses
        try
            x_sol = fsolve(ode_fun, x0, optimoptions('fsolve', 'Display', 'none', ...
                'TolFun', 1e-14, 'TolX', 1e-14, 'MaxIterations', 10000));
            if ~isnan(x_sol) && ~isinf(x_sol) && x_sol > 0 && x_sol < 1
                gamma_x = beta - (beta - alpha) * x_sol;
                if gamma_x > 1e-8
                    H_sol = (B * x_sol / (rho + x_sol) - d) / gamma_x;
                    if ~isnan(H_sol) && ~isinf(H_sol)
                        residual = abs(ode_fun(x_sol));
                        if residual < 1e-10 && (isempty(x_sols) || min(abs(x_sols - x_sol)) > 1e-6)
                            x_sols = [x_sols; x_sol];
                            H_sols = [H_sols; H_sol];
                            fprintf('Root candidate: x=%.12f, H=%.12f, gamma_x=%.12f, residual=%.2e\n', ...
                                    x_sol, H_sol, gamma_x, residual);
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
    for i = 1:min(length(x_sols), 4)
        eq_points(2+i, :) = [x_sols(i), H_sols(i)];
        x_eq = x_sols(i);
        H_eq = H_sols(i);
        gamma_x = beta - (beta - alpha) * x_eq;
        % Jacobian: d(dH/dt)/dH, d(dH/dt)/dx, d(dx/dt)/dH, d(dx/dt)/dx
        J = [(B * x_eq / (rho + x_eq) - d - 2 * gamma_x * H_eq), ...
             (B * H_eq * rho) / (rho + x_eq)^2 + H_eq^2 * (beta - alpha); ...
             -c * x_eq / (rho + x_eq), ...
             1 - 2 * x_eq - c * H_eq * rho / (rho + x_eq)^2];
        fprintf('Jacobian at (x=%.12f, H=%.12f):\n', x_eq, H_eq);
        fprintf('[%.12f, %.12f]\n[%.12f, %.12f]\n', J(1,1), J(1,2), J(2,1), J(2,2));
        eigvals = eig(J);
        fprintf('Eigenvalues real parts: [%.12f, %.12f], imaginary parts: [%.12f, %.12f]\n', ...
                real(eigvals(1)), real(eigvals(2)), imag(eigvals(1)), imag(eigvals(2)));
        is_stable(2+i) = all(real(eigvals) < -1e-8);
        eq_types{2+i} = classify_equilibrium(eigvals, is_stable(2+i));
        fprintf('Equilibrium %d (x=%.12f, H=%.12f): Eigenvalues = [%.12f, %.12f], Stable=%d, Type=%s\n', ...
                2+i, x_eq, H_eq, eigvals(1), eigvals(2), is_stable(2+i), eq_types{2+i});
    end
end

function type = classify_equilibrium(eigvals, is_stable)
    if any(abs(imag(eigvals)) > 1e-8)
        if is_stable
            type = 'Stable Spiral';
        else
            type = 'Unstable Spiral';
        end
    else
        if eigvals(1) * eigvals(2) < -1e-8
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