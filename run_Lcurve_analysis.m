function [lambda_optimal] = run_l_curve_analysis(exp_tmv_V, thetas, nu, N, j_val, pole_type, initial_guesses, bounds, optim_options)
%% 【鲁棒版】运行L曲线分析，寻找最佳的lambda值
    fprintf('--- 开始L曲线分析以寻找最佳lambda ---\n');
    
    lambdas = logspace(-12, -2, 40);
    residual_norms = zeros(size(lambdas));
    solution_norms = zeros(size(lambdas));
    
    h = waitbar(0, '正在进行L曲线分析...');
    for i = 1:length(lambdas)
        lambda_current = lambdas(i);
        
        best_err_fit = inf;
        best_params = [];
        for guess_idx = 1:size(initial_guesses, 1)
            handle = @(p) objective_function_tmv_regularized(p, exp_tmv_V, thetas, nu, N, j_val, pole_type, lambda_current);
            [params, ~] = fmincon(handle, initial_guesses(guess_idx, :), [], [], [], [], bounds.lower, bounds.upper, [], optim_options);
            [err_fit, ~] = objective_function_tmv_regularized(params, exp_tmv_V, thetas, nu, N, j_val, pole_type, 0);
            if err_fit < best_err_fit
                best_err_fit = err_fit;
                best_params = params;
            end
        end
        residual_norms(i) = best_err_fit;
        solution_norms(i) = best_params(1)^2;
        waitbar(i/length(lambdas), h, sprintf('测试 Lambda: %.2e', lambda_current));
    end
    close(h);
    
    figure('Name', '图1：L曲线分析', 'NumberTitle', 'off');
    loglog(residual_norms, solution_norms, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    xlabel('残差范数 (拟合误差)');
    ylabel('解的范数 (惩罚项)');
    title(['L曲线分析 (' pole_type ' 极)']);
    grid on;
    
    x = log(residual_norms); y = log(solution_norms);
    dx = gradient(x); dy = gradient(y);
    ddx = gradient(dx); ddy = gradient(dy);
    curvature = abs(dx.*ddy - dy.*ddx) ./ (dx.^2 + dy.^2).^(3/2);
    search_range = floor(length(curvature)*0.2) : floor(length(curvature)*0.8);
    [~, max_curv_idx_local] = max(curvature(search_range));
    max_curv_idx = search_range(max_curv_idx_local);
    lambda_optimal = lambdas(max_curv_idx);
    
    fprintf('基于最大曲率法自动选择的最佳Lambda: %.4e\n', lambda_optimal);
    hold on;
    loglog(residual_norms(max_curv_idx), solution_norms(max_curv_idx), 'r*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', sprintf('自动选择 λ = %.2e', lambda_optimal));
    
    indices_to_label = round(linspace(1, length(lambdas), 8));
    for i = indices_to_label
        text(residual_norms(i), solution_norms(i)*1.5, sprintf('λ=%.1e', lambdas(i)), 'FontSize', 9, 'Color', 'k');
    end
    
    legend('show', 'Location', 'best');
    set(gca, 'FontSize', 12);
    fprintf('请检查L曲线图。如果自动选择的拐角不理想，请在主程序中手动覆盖lambda值。\n');
end