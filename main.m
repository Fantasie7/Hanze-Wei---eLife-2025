%% 【最终版主程序 - 从文件读取】
% 该程序从指定的Excel文件读取TMV数据，并对每一行数据（时间点）
% 进行完整的正则化分析与对比可视化。

clear; clc; close all;

%% 1. 全局参数和文件路径设置
% --- 用户需要修改的部分 ---
file_path = ''; % <--- 请确认您的Excel文件路径
sheets_to_process = {'5'}; % <--- 请指定要处理的工作表，例如 {'5', '12'}
% --- 参数设置 ---
nu = 0.5;        % 电导率比 (sigma_i / sigma_e)
N = 59;           % 勒让德多项式展开的最高阶数
j_val = 0,1,2;       % Gm(θ) 形状参数

% 为优化器提供多组初始猜测值 [g0, theta_c_deg]
initial_guesses_pole = [1, 20; 10, 45; 50, 70; 100, 85];
% 为参数 [g0, theta_c] 定义边界
bounds.lower = [1e-4, 0];
bounds.upper = [100, 90];
optim_options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');

% --- 角度映射与分组 ---
all_theta_mapping = [60, 30, 0, 30, 60, 90, 120, 150, 180, 150, 120, 90];
anode_indices = 1:6;
cathode_indices = 7:12;
anode_thetas_deg = all_theta_mapping(anode_indices);
cathode_thetas_deg = all_theta_mapping(cathode_indices);

%% 2. 循环处理每个工作表和每一行数据
for i = 1:length(sheets_to_process)
    sheet_name = sheets_to_process{i};
    fprintf('================================================\n');
    fprintf('开始处理工作表: %s kV/cm\n', sheet_name);
    
    data = readtable(file_path, 'Sheet', sheet_name, 'VariableNamingRule', 'preserve');
    num_rows = height(data);

    for t_idx = 1:num_rows
        time_val_ns = data{t_idx, 1};
        exp_tmv_mV_all = data{t_idx, 2:2:25};
        
        if any(isnan(exp_tmv_mV_all))
            fprintf('  - 跳过 t = %.1f ns，数据不完整。\n', time_val_ns);
            continue;
        end
        
        fprintf('\n--- 正在分析时间点 t = %.1f ns ---\n', time_val_ns);
        
        % 分割阳极和阴极的实验数据
        exp_tmv_V_anode = exp_tmv_mV_all(anode_indices) / 1000;
        exp_tmv_V_cathode = exp_tmv_mV_all(cathode_indices) / 1000;
        
        % =================== A. 分析阳极 (Anode) ===================
        fprintf('\n--- (1/2) 正在分析阳极 (Anode) ---\n');
        pole_type = 'anode';
        
        % A1. L曲线分析确定最佳lambda
        lambda_opt_anode = run_l_curve_analysis(exp_tmv_V_anode, anode_thetas_deg, nu, N, j_val, pole_type, initial_guesses_pole, bounds, optim_options);
        
        % --- 手动调整通道 ---
        % lambda_manual = 1e-6; % <--- 如果需要，在此处手动设置lambda
        % lambda_opt_anode = lambda_manual;
        % fprintf('注意：已手动覆盖阳极Lambda值为 %.2e\n', lambda_manual);
        
        fprintf('\n将使用 lambda_opt = %.4e 对阳极数据进行最终拟合。\n', lambda_opt_anode);

        % A2. 进行有/无正则化的拟合
        [params_nr, ~] = perform_fitting('non-regularized', exp_tmv_V_anode, anode_thetas_deg, nu, N, j_val, pole_type, 0, initial_guesses_pole, bounds, optim_options);
        [params_r, ~] = perform_fitting('regularized', exp_tmv_V_anode, anode_thetas_deg, nu, N, j_val, pole_type, lambda_opt_anode, initial_guesses_pole, bounds, optim_options);
        
        % A3. 可视化对比结果
        visualize_comparison_results(exp_tmv_V_anode, params_nr, params_r, anode_thetas_deg, nu, N, j_val, pole_type, time_val_ns, sheet_name);

        % =================== B. 分析阴极 (Cathode) ===================
        fprintf('\n--- (2/2) 正在分析阴极 (Cathode) ---\n');
        pole_type = 'cathode';

        % B1. L曲线分析确定最佳lambda
        lambda_opt_cathode = run_l_curve_analysis(exp_tmv_V_cathode, cathode_thetas_deg, nu, N, j_val, pole_type, initial_guesses_pole, bounds, optim_options);

        % B2. 进行有/无正则化的拟合
        [params_nr, ~] = perform_fitting('non-regularized', exp_tmv_V_cathode, cathode_thetas_deg, nu, N, j_val, pole_type, 0, initial_guesses_pole, bounds, optim_options);
        [params_r, ~] = perform_fitting('regularized', exp_tmv_V_cathode, cathode_thetas_deg, nu, N, j_val, pole_type, lambda_opt_cathode, initial_guesses_pole, bounds, optim_options);

        % B3. 可视化对比结果
        visualize_comparison_results(exp_tmv_V_cathode, params_nr, params_r, cathode_thetas_deg, nu, N, j_val, pole_type, time_val_ns, sheet_name);
        
        fprintf('\n--- 时间点 t = %.1f ns 分析完成 ---\n', time_val_ns);
        if t_idx < num_rows
            pause; % 暂停，按任意键继续处理下一行数据
        end
    end
end
fprintf('\n================================================\n');
fprintf('所有指定工作表和数据行已处理完毕。\n');


%% --- 辅助函数 ---
function [best_params, best_err] = perform_fitting(type, exp_tmv, thetas, nu, N, j, pole, lambda, guesses, bounds, options)
    fprintf('--- 正在进行【%s】拟合... ---\n', type);
    best_err = inf; best_params = [];
    for i = 1:size(guesses, 1)
        if strcmp(type, 'regularized')
            handle = @(p) objective_function_tmv_regularized(p, exp_tmv, thetas, nu, N, j, pole, lambda);
            [params, ~] = fmincon(handle, guesses(i,:), [], [], [], [], bounds.lower, bounds.upper, [], options);
            [err_fit, ~] = objective_function_tmv_regularized(params, exp_tmv, thetas, nu, N, j, pole, 0);
             if err_fit < best_err, best_err = err_fit; best_params = params; end
        else % non-regularized
            handle = @(p) objective_function_tmv(p, exp_tmv, thetas, nu, N, j, pole);
            [params, err] = fmincon(handle, guesses(i,:), [], [], [], [], bounds.lower, bounds.upper, [], options);
            if err < best_err, best_err = err; best_params = params; end
        end
    end
    fprintf('完成！结果: g0 = %.2f, tc = %.1f\n', best_params(1), best_params(2));
end

function visualize_comparison_results(exp_tmv_noisy, p_nr, p_r, thetas_deg, nu, N, j, pole, time, sheet)
    % 与演示版不同，这里没有“真实值”，只对比有/无正则化的结果
    placeholder = zeros(size(thetas_deg));
    [~, tmv_fit_nr] = objective_function_tmv(p_nr, placeholder, thetas_deg, nu, N, j, pole);
    [~, tmv_fit_r] = objective_function_tmv(p_r, placeholder, thetas_deg, nu, N, j, pole);
    
    fig_title_prefix = sprintf('%s kV/cm, t=%.1fns, %s 极', sheet, time, pole);

    % 图2：TMV拟合对比
    figure('Name', sprintf('TMV对比 @ t=%.1fns (%s)', time, pole), 'NumberTitle', 'off');
    hold on;
    plot(thetas_deg, abs(exp_tmv_noisy), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', '实验数据');
    plot(thetas_deg, abs(tmv_fit_nr), 'r--s', 'LineWidth', 2, 'DisplayName', '无正则化拟合');
    plot(thetas_deg, abs(tmv_fit_r), 'b-^', 'LineWidth', 2, 'DisplayName', '有正则化拟合');
    hold off; grid on; xlabel('角度 (度)'); ylabel('跨膜电压 |TMV| (V)');
    title({fig_title_prefix, 'TMV数据拟合效果对比'});
    legend('show', 'Location', 'best'); set(gca, 'FontSize', 12);

    % 图3：Gm重构对比
    theta_plot_deg = linspace(min(thetas_deg), max(thetas_deg), 200);
    th_rad = deg2rad(theta_plot_deg);
    
    gm_nr = p_nr(1) * ((abs(cos(th_rad)) - cosd(p_nr(2))).^j) ./ ((1 - cosd(p_nr(2))).^j) .* (abs(cos(th_rad)) > cosd(p_nr(2)));
    gm_r = p_r(1) * ((abs(cos(th_rad)) - cosd(p_r(2))).^j) ./ ((1 - cosd(p_r(2))).^j) .* (abs(cos(th_rad)) > cosd(p_r(2)));

    figure('Name', sprintf('Gm对比 @ t=%.1fns (%s)', time, pole), 'NumberTitle', 'off');
    hold on;
    plot(theta_plot_deg, gm_nr, 'r--', 'LineWidth', 2, 'DisplayName', '无正则化重构 (不稳定)');
    plot(theta_plot_deg, gm_r, 'b-', 'LineWidth', 2, 'DisplayName', '有正则化重构 (稳定)');
    hold off; grid on; xlabel('角度 (度)'); ylabel('膜电导 Gm (S/m^2)');
    title({fig_title_prefix, '最终Gm(θ)重构结果对比'});
    legend('show', 'Location', 'northeast'); set(gca, 'FontSize', 12);
end