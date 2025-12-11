%% 主程序：基于H1991的TMV无量纲化拟合（E=5/12/50/90 kV/cm）
% 关键修改：1. TMV无量纲化（Δψ=ΔΨ/(aE0)）；2. 绑定sheet与场强E0；3. 修正初始值与文献匹配
clear; clc; close all;

%% 1. 全局参数设置（严格匹配文献定义）
% 基础物理参数（文献附录+用户指定）
a_cm = 5e-4; % 细胞半径：5μm = 5×10^-4 cm（用户指定）
nu = 0.5; % 胞内/胞外电导率比（κ_i/κ_e）
N = 59;  % 勒让德多项式截断阶数
j_values_to_test = [0,1,2]; % 完整测试文献中的j=0/1/2

% 【核心】绑定sheet名称与对应的场强E0（单位：V/cm，kV/cm→V/cm换算）
sheet_E0_map = containers.Map({'5', '12', '50', '90'},[5000, 12000, 50000, 90000]);

% 无量纲电导g0初始猜测
initial_guesses_pole =... 
[
    0.1, 30;   % 
    2.0, 60;   % 
    4.0, 80    % 
];

% 参数边界（基于文献范围）
lower_bounds = [0.01, 0];   % g0下限（文献最小值），θc下限0°
upper_bounds = [200, 90];    % g0上限（文献最大值），θc上限90°
optim_options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', 'MaxFunEvals', 5000);

% 角度映射（文献Fig.2-4：0°=阳极，180°=阴极，90°=赤道）
all_theta_mapping = [60, 30, 0, 30, 60, 90, 120, 150, 180, 150, 120, 90];
anode_indices = 1:6;   % 阳极区域：θ=0°/30°/60°/90°（含赤道，非穿孔区）
cathode_indices = 7:12;% 阴极区域：θ=90°/120°/150°/180°（含赤道，非穿孔区）
anode_thetas = all_theta_mapping(anode_indices);
cathode_thetas = all_theta_mapping(cathode_indices);

% 文件路径
file_path = '';
sheets_to_process = {'5', '12', '50', '90'}; % 待处理的场强sheet

%% 2. 循环处理每个Sheet（按场强匹配E0）
output_folder = 'Fitting_Results_Dimensionless'; 
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

for i = 1:length(sheets_to_process)
    sheet_name = sheets_to_process{i};
    E0_Vcm = sheet_E0_map(sheet_name); % 读取当前sheet对应的场强E0（V/cm）
    fprintf('================================================\n');
    fprintf('处理 Sheet: %s kV/cm | E0=%.0f V/cm\n', sheet_name, E0_Vcm);
    
    % 初始化输出文件
    log_filename = fullfile(output_folder, sprintf('Fitting_Log_%skV.xlsx', sheet_name));
    cond_filename = fullfile(output_folder, sprintf('Conductivity_Results_%skV.xlsx', sheet_name));
    param_summary_filename = fullfile(output_folder, sprintf('Parameter_Summary_%skV.xlsx', sheet_name));
    if isfile(log_filename), delete(log_filename); end
    if isfile(param_summary_filename), delete(param_summary_filename); end

    % 读取实验数据（TMV单位：mV）
    data = readtable(file_path, 'Sheet', sheet_name, 'VariableNamingRule', 'preserve');
    num_rows = height(data);
    
    summary_table = data(:, 1);
    summary_table.best_j = NaN(num_rows, 1); % 修正：阴阳极共用j（文献对称性）
    summary_table.g0 = NaN(num_rows, 1);     % 修正：阴阳极共用g0（文献对称性）
    summary_table.theta_c = NaN(num_rows, 1); % 修正：阴阳极共用θc（文献对称性）
    conductivity_data = NaN(num_rows, 12);

    % 逐行拟合（修正：阴阳极不独立，按文献对称拟合）
    fprintf('开始无量纲化拟合...\n');
    for t_idx = 1:num_rows
        time_val_ns = data{t_idx, 1};
        exp_tmv_mV_all = data{t_idx, 2:2:25}; % 实验TMV（mV）
        
        if any(isnan(exp_tmv_mV_all))
            fprintf('  - 跳过 t=%.1f ns：数据不完整\n', time_val_ns);
            continue;
        end
        
        % 步骤1：实际TMV单位转换（mV→V）
        exp_tmv_V_all = exp_tmv_mV_all / 1000;
        
        % 步骤2：TMV无量纲化（核心！按文献Δψ=ΔΨ/(a*E0)）
        dimensionless_factor = a_cm * E0_Vcm; % a*E0（单位：cm*V/cm=V）
        exp_tmv_psi_all = exp_tmv_V_all / dimensionless_factor; % 无量纲Δψ
        
        % 分割阴阳极无量纲TMV（用于拟合输入）
        exp_tmv_psi_anode = exp_tmv_psi_all(anode_indices);
        exp_tmv_psi_cathode = exp_tmv_psi_all(cathode_indices);
        % 合并阴阳极数据（文献对称性要求：拟合统一参数，不独立拆分）
        exp_tmv_psi_combined = [exp_tmv_psi_anode, exp_tmv_psi_cathode];
        combined_thetas = [anode_thetas, cathode_thetas];

        % 步骤3：拟合（修正：阴阳极共用参数，符合文献对称性）
        best_error = inf;
        best_params = [];
        best_j = NaN;
        for j_val = j_values_to_test
            for guess_idx = 1:size(initial_guesses_pole, 1)
                % 目标函数输入：无量纲TMV（exp_tmv_psi_combined）
                handle = @(p) objective_function_tmv(p, exp_tmv_psi_combined, combined_thetas, nu, N, j_val);
                [params, err] = fmincon(handle, initial_guesses_pole(guess_idx, :), [], [], [], [], lower_bounds, upper_bounds, [], optim_options);
                if err < best_error
                    best_error = err;
                    best_params = params;
                    best_j = j_val;
                end
            end
        end

        fprintf('  - t=%.1f ns | j=%d, g0=%.2f, θc=%.1f° | 无量纲拟合误差=%.4f\n', ...
            time_val_ns, best_j, best_params(1), best_params(2), best_error);

        % 步骤4：结果记录
        summary_table.best_j(t_idx) = best_j;
        summary_table.g0(t_idx) = best_params(1);
        summary_table.theta_c(t_idx) = best_params(2);

        % 计算理论无量纲TMV与实际电导
        [~, theory_tmv_psi] = objective_function_tmv(best_params, exp_tmv_psi_combined, combined_thetas, nu, N, best_j);
        % 计算实际膜电导G（文献Eq.：G=(g*κ_i)/a，κ_i=nu*κ_e，取κ_e=5 S/m=0.05 S/cm）
        kappa_e_Scm = 0.05; % 海水电导率（文献常用值）
        kappa_i_Scm = nu * kappa_e_Scm;
        g0 = best_params(1); % 无量纲电导
        theta_c_rad = deg2rad(best_params(2));
        j_val = best_j;
        g_theta = zeros(1, 12); % 实际电导G（单位：S/cm²）
        for k = 1:12
            th_rad = deg2rad(all_theta_mapping(k));
            % 穿孔区域判断（文献Eq.A16：θ≤θc 或 θ≥180°-θc）
            is_porated = (th_rad <= theta_c_rad) || (th_rad >= pi - theta_c_rad);
            if is_porated
                if j_val == 0
                    g_dimless = g0;
                else
                    numerator = (abs(cos(th_rad)) - cos(theta_c_rad))^j_val;
                    denominator = (1 - cos(theta_c_rad))^j_val;
                    g_dimless = g0 * numerator / denominator;
                end
                G_Scm2 = (g_dimless * kappa_i_Scm) / a_cm; % 实际电导（S/cm²）
                g_theta(k) = G_Scm2 * 100; % 转换为S/m²（1 S/cm²=100 S/m²）
            else
                g_theta(k) = 0; % 非穿孔区电导为0
            end
        end
        conductivity_data(t_idx, :) = g_theta;

        % 日志记录（含无量纲与实际值对照）
        log_table = table(combined_thetas', exp_tmv_psi_combined', theory_tmv_psi', ...
            'VariableNames', {'Angle_deg', 'Exp_Dimensionless_Δψ', 'Theory_Dimensionless_Δψ'});
        % 转换为实际TMV（V）用于直观查看
        log_table.Exp_Actual_TMV_V = log_table.Exp_Dimensionless_Δψ * dimensionless_factor;
        log_table.Theory_Actual_TMV_V = log_table.Theory_Dimensionless_Δψ * dimensionless_factor;
        sheet_name_log = sprintf('t_%.1fns', time_val_ns);
        sheet_name_log = strrep(sheet_name_log, '-', 'neg_');
        writetable(log_table, log_filename, 'Sheet', sheet_name_log);
    end

    %% 3. 结果可视化（修正：标注无量纲/实际单位）
    fprintf('生成结果图像...\n');
    time_vector = summary_table{:, 1};
    num_time_points = length(time_vector);
    num_angle_points = 12;

    % 图1：电导率 vs 角度（实际单位：S/m²）
    figure('Name', sprintf('Conductivity vs Angle (%s kV/cm)', sheet_name), 'Visible', 'off');
    hold on; colors = jet(num_time_points);
    legend_entries = cell(num_time_points, 1);
    for t_idx = 1:num_time_points
        plot(1:num_angle_points, conductivity_data(t_idx, :), '-o', 'LineWidth', 1.5, 'Color', colors(t_idx,:));
        legend_entries{t_idx} = sprintf('t=%.1f ns', time_vector(t_idx));
    end
    title(sprintf('膜电导率 vs 角度位置（%s kV/cm，单位：S/m²）', sheet_name));
    xlabel('角度位置索引（1-12，对应θ=60°~90°~180°）');
    ylabel('实际膜电导率（S/m²）');
    xticks(1:num_angle_points); grid on;
    legend(legend_entries, 'Location', 'northeastoutside');
    plot1_filename = fullfile(output_folder, sprintf('Cond_vs_Angle_%skV.png', sheet_name));
    saveas(gcf, plot1_filename);
    fprintf('  - 保存图像：%s\n', plot1_filename);

    % 图2：无量纲Δψ vs 角度（文献标准输出）
    figure('Name', sprintf('Dimensionless Δψ vs Angle (%s kV/cm)', sheet_name), 'Visible', 'off');
    hold on; colors = jet(num_time_points);
    for t_idx = 1:num_time_points
        exp_tmv_V_all = data{t_idx, 2:2:25}/1000;
        exp_tmv_psi_all = exp_tmv_V_all / (a_cm * E0_Vcm);
        plot(1:num_angle_points, exp_tmv_psi_all, 'o-', 'LineWidth', 1.5, 'Color', colors(t_idx,:));
    end
    title(sprintf('无量纲跨膜电位Δψ vs 角度位置（%s kV/cm）', sheet_name));
    xlabel('角度位置索引（1-12）');
    ylabel('无量纲跨膜电位Δψ（无单位）');
    xticks(1:num_angle_points); grid on;
    legend(legend_entries, 'Location', 'northeastoutside');
    plot2_filename = fullfile(output_folder, sprintf('DeltaPsi_vs_Angle_%skV.png', sheet_name));
    saveas(gcf, plot2_filename);
    fprintf('  - 保存图像：%s\n', plot2_filename);

    close all;

    % 保存结果文件
    writetable(summary_table, param_summary_filename);
    conductivity_results_table = data(:, 1);
    for k = 1:12
        conductivity_results_table.(sprintf('Loc_%d_Conductivity_Sm2', k)) = conductivity_data(:, k);
    end
    writetable(conductivity_results_table, cond_filename);
    fprintf('  - 保存参数文件：%s\n', param_summary_filename);
    fprintf('  - 保存电导率文件：%s\n', cond_filename);
end
fprintf('================================================\n');
fprintf('所有场强拟合完成！结果路径：%s\n', output_folder);


%% 目标函数（修正：输入输出均为无量纲Δψ，符合文献）
function [error, theory_tmv_psi] = objective_function_tmv(params, exp_tmv_psi, thetas_deg, nu, N, j_val)
% 输入：exp_tmv_psi - 无量纲实验TMV（Δψ）
% 输出：theory_tmv_psi - 无量纲理论TMV（Δψ），error - 无量纲拟合误差

    % 解包参数（g0=无量纲电导，θc=穿孔临界角度）
    g0 = params(1);
    theta_c_rad = deg2rad(params(2));
    j = j_val;

    % 1. 膜电导函数（文献Eq.A16：无量纲电导g(θ)）
    g_func = @(th_rad) begin
        is_porated = (th_rad <= theta_c_rad) || (th_rad >= pi - theta_c_rad);
        if ~is_porated
            0;
        else
            if j == 0
                g0;
            else
                numerator = (abs(cos(th_rad)) - cos(theta_c_rad))^j;
                denominator = (1 - cos(theta_c_rad))^j;
                g0 * numerator / denominator;
            end
        end
    end;
    g_func = vectorize(g_func); % 向量化函数

    % 2. 计算g_mn矩阵（文献Eq.A14：无量纲电导的积分）
    odd_indices = 1:2:N; % 仅奇数阶（文献对称性要求）
    num_odd = length(odd_indices);
    g_mn = zeros(num_odd, num_odd);
    theta_int = linspace(0, pi, 200); % 积分步长（200点，保证精度）
    x_int = cos(theta_int);
    for m_idx = 1:num_odd
        m = odd_indices(m_idx);
        Pm = legendre_poly(m, x_int);
        for n_idx = 1:num_odd
            n = odd_indices(n_idx);
            Pn = legendre_poly(n, x_int);
            g_vals = g_func(theta_int);
            integrand = Pm .* g_vals .* Pn;
            g_mn(m_idx, n_idx) = trapz(x_int, integrand); % 数值积分
        end
    end

    % 3. 构建线性方程组（文献Eq.A13：A*c = b）
    A = zeros(num_odd, num_odd);
    b = zeros(num_odd, 1);
    g_m1 = g_mn(:, 1); % g_m1 = ∫Pm*g*P1 dx（P1为1阶勒让德多项式）
    for m_idx = 1:num_odd
        m = odd_indices(m_idx);
        % 右侧b：3/2 * g_m1（文献Eq.A13）
        b(m_idx) = (3/2) * g_m1(m_idx);
        % 左侧A：对角线项+非对角线项（修正原程序重复计算错误）
        for n_idx = 1:num_odd
            n = odd_indices(n_idx);
            coeff = (nu/(n+1) + 1/n) * g_mn(m_idx, n_idx);
            if m_idx == n_idx
                A(m_idx, n_idx) = 2/(2*m + 1) + coeff; % 对角线项：2/(2m+1) + coeff
            else
                A(m_idx, n_idx) = coeff; % 非对角线项：仅coeff
            end
        end
    end

    % 4. 求解勒让德系数c（文献Eq.A9-A10）
    c_coeffs = A \ b; % 线性方程组求解

    % 5. 计算无量纲理论TMV（文献Eq.：Δψ=3/2 P1 - Σc_n(ν/(n+1)+1/n)Pn）
    x_exp = cosd(thetas_deg);
    theory_tmv_psi = (3/2) * legendre_poly(1, x_exp); % 基础项（3/2 P1）
    for idx = 2:num_odd
        n = odd_indices(idx);
        Pn = legendre_poly(n, x_exp);
        theory_tmv_psi = theory_tmv_psi - c_coeffs(idx) * (nu/(n+1) + 1/n) * Pn;
    end

    % 6. 无量纲误差计算（文献：基于Δψ的均方根误差）
    error = sqrt(mean((theory_tmv_psi - exp_tmv_psi).^2));
    if ~isreal(error) || isnan(error), error = 1e9; end
end


%% 勒让德多项式计算（文献标准公式，修正原程序无错误）
function Pn = legendre_poly(n, x)
% 文献附录推荐：n阶勒让德多项式（Rodrigues公式）
    if n < 0 || ~isinteger(n)
        error('勒让德多项式阶数n必须为非负整数（文献附录要求）');
    end
    Pn = 0;
    k_max = floor(n/2);
    for k = 0:k_max
        numerator = (-1)^k * factorial(2*n - 2*k);
        denominator = 2^n * factorial(k) * factorial(n - k) * factorial(n - 2*k);
        term = numerator / denominator * (x .^ (n - 2*k));
        Pn = Pn + term;
    end
end


% %% 主程序：通过独立拟合阴阳极来计算膜电导
% % 新功能描述:
% %   1. 【核心】对阳极(6点)和阴极(6点)分别进行独立的参数拟合。
% %   2. 使用多组初始值进行尝试，提高结果的稳健性。
% %   3. 保持对 j=0, 1, 2 的对比，为阴阳极分别寻找最佳j值。
% %   4. 输出格式保持不变，但结果基于更精确的拟合。
% clear; clc; close all;
% 
% %% 1. 全局参数设置
% nu = 0.5; % nu = sigma_i/sigma_e ; sigma_i = 0.5 ; sigma_e = 1.0
% N = 59;  %
% j_values_to_test =  [0,1,2]
% 
% % --- 多组初始猜测值 [g0, theta_c_deg] ---
% % 这个列表将同时用于阳极和阴极的拟合
% initial_guesses_pole =... 
% [
% ];
% 
% % 为参数 [g0, theta_c] 定义边界
% optim_options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', 'MaxFunEvals', 1000);
% 
% % --- 角度映射与分组 ---
% all_theta_mapping = [60, 30, 0, 30, 60, 90, 120, 150, 180, 150, 120, 90];
% anode_indices = 1:6;
% cathode_indices = 7:12;
% anode_thetas = all_theta_mapping(anode_indices);
% cathode_thetas = all_theta_mapping(cathode_indices);
% 
% % --- 文件路径 ---
% file_path = '';
% sheets_to_process = {'5', '12', '50', '90'}; % , '100'
% 
% %% 2. 循环处理每个Sheet
% % 定义输出文件夹的名称
% output_folder = 'Fitting_Results'; 
% %检查文件夹是否存在，如果不存在则创建它
% if ~exist(output_folder, 'dir')
%    mkdir(output_folder);
%    fprintf('已创建输出文件夹: %s\n', output_folder);
% end
% 
% for i = 1:length(sheets_to_process)
%     sheet_name = sheets_to_process{i};
%     fprintf('================================================\n');
%     fprintf('开始处理 Sheet: %s kV/cm\n', sheet_name);
% 
%     % --- 初始化 ---
%     log_filename = fullfile(output_folder, sprintf('Fitting_Log_%skV.xlsx', sheet_name));
%     cond_filename = fullfile(output_folder, sprintf('Conductivity_Results_%skV.xlsx', sheet_name));
%     param_summary_filename = fullfile(output_folder, sprintf('Parameter_Summary_%skV.xlsx', sheet_name));
%     if isfile(log_filename), delete(log_filename); end
%     if isfile(param_summary_filename), delete(param_summary_filename); end
% 
%     data = readtable(file_path, 'Sheet', sheet_name, 'VariableNamingRule', 'preserve');
%     num_rows = height(data);
% 
%     summary_table = data(:, 1);
%     % 扩展摘要表以存储6个参数
%     summary_table.best_j_anode = NaN(num_rows, 1);
%     summary_table.g0_anode = NaN(num_rows, 1);
%     summary_table.theta_c_anode = NaN(num_rows, 1);
%     summary_table.best_j_cathode = NaN(num_rows, 1);
%     summary_table.g0_cathode = NaN(num_rows, 1);
%     summary_table.theta_c_cathode = NaN(num_rows, 1);
%     conductivity_data = NaN(num_rows, 12);
% 
%     % --- 对sheet进行逐行拟合 ---
%     fprintf('开始对Sheet "%s" 进行逐行、独立的阴阳极拟合...\n', sheet_name);
%     for t_idx = 1:num_rows
%         time_val_ns = data{t_idx, 1};
%         exp_tmv_mV_all = data{t_idx, 2:2:25};
% 
%         if any(isnan(exp_tmv_mV_all))
%             fprintf('  - 跳过 t = %.1f ns，数据不完整。\n', time_val_ns);
%             continue;
%         end
% 
%         % 分割阴极和阳极跨膜电压数据
%         exp_tmv_V_anode = exp_tmv_mV_all(anode_indices) / 1000;     %这里没有进行无量纲化处理，请按不同E 5 12 50 90 kV/cm处理对应的sheet
%         exp_tmv_V_cathode = exp_tmv_mV_all(cathode_indices) / 1000;
% 
%         % --- 【A. 独立拟合阳极】---
%         best_anode_error = inf;
%         best_anode_params = [];
%         best_anode_j = NaN;
%         for j_val = j_values_to_test
%             for guess_idx = 1:size(initial_guesses_pole, 1)
%                 handle = @(p) objective_function_tmv(p, exp_tmv_V_anode, anode_thetas, nu, N, j_val, 'anode');
%                 % --- 使用 fmincon 进行优化 ---
%                 % fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options)
%                 % 只使用 x0, lb, ub，其他线性约束设为空[]
%                 [params, err] = fmincon(handle, initial_guesses_pole(guess_idx, :), [], [], [], [], lower_bounds, upper_bounds, [], optim_options);
%                 if err < best_anode_error
%                     best_anode_error = err;
%                     best_anode_params = params;
%                     best_anode_j = j_val;
%                 end
%             end
%         end
% 
%         % --- 【B. 独立拟合阴极】---
%         best_cathode_error = inf;
%         best_cathode_params = [];
%         best_cathode_j = NaN;
%         for j_val = j_values_to_test
%             for guess_idx = 1:size(initial_guesses_pole, 1)
%                 handle = @(p) objective_function_tmv(p, exp_tmv_V_cathode, cathode_thetas, nu, N, j_val, 'cathode');
%                 % --- 使用 fmincon 进行优化 ---
%                 [params, err] = fmincon(handle, initial_guesses_pole(guess_idx, :), [], [], [], [], lower_bounds, upper_bounds, [], optim_options);
%                 if err < best_cathode_error
%                     best_cathode_error = err;
%                     best_cathode_params = params;
%                     best_cathode_j = j_val;
%                 end
%             end
%         end
% 
%         fprintf('  - 完成 t=%.1f ns | 阳极(j=%d, g0=%.2f, tc=%.1f) | 阴极(j=%d, g0=%.2f, tc=%.1f)\n', ...
%             time_val_ns, best_anode_j, best_anode_params(1), best_anode_params(2), ...
%             best_cathode_j, best_cathode_params(1), best_cathode_params(2));
% 
%         % --- 【C. 结果汇总与记录】---
%         % 存储参数摘要
%         summary_table.best_j_anode(t_idx) = best_anode_j;
%         summary_table.g0_anode(t_idx) = best_anode_params(1);
%         summary_table.theta_c_anode(t_idx) = best_anode_params(2);
%         summary_table.best_j_cathode(t_idx) = best_cathode_j;
%         summary_table.g0_cathode(t_idx) = best_cathode_params(1);
%         summary_table.theta_c_cathode(t_idx) = best_cathode_params(2);
% 
%         % 分别计算理论TMV并合并，用于日志记录
%         [~, theory_tmv_anode] = objective_function_tmv(best_anode_params, exp_tmv_V_anode, anode_thetas, nu, N, best_anode_j, 'anode');
%         [~, theory_tmv_cathode] = objective_function_tmv(best_cathode_params, exp_tmv_V_cathode, cathode_thetas, nu, N, best_cathode_j, 'cathode');
%         best_theory_tmv_all = [theory_tmv_anode, theory_tmv_cathode];
% 
%         log_table = table(all_theta_mapping', exp_tmv_mV_all'/1000, abs(best_theory_tmv_all)', 'VariableNames', {'Angle_deg', 'Experimental_TMV_V', 'Best_Fit_Theoretical_TMV_V'});
%         log_table.Difference = abs(log_table.Experimental_TMV_V - log_table.Best_Fit_Theoretical_TMV_V);
%         sheet_name_log = sprintf('t_%.1fns', time_val_ns);
%         sheet_name_log = strrep(sheet_name_log, '-', 'neg_');
%         writetable(log_table, log_filename, 'Sheet', sheet_name_log);
% 
%         % 计算最终电导率
%         % 此处使用独立拟合出的最佳参数，为所有12个位点计算最终的膜电导率
%         g_theta = zeros(1, 12); % 初始化电导率向量，默认为0（未穿孔区域）
% 
%         g0_a = best_anode_params(1);
%         tc_a_rad = deg2rad(best_anode_params(2)); 
%         j_a = best_anode_j;
% 
%         g0_c = best_cathode_params(1);
%         tc_c_rad = deg2rad(best_cathode_params(2)); 
%         j_c = best_cathode_j;
% 
%         % 遍历所有12个位点，计算当前位点电导率
%         for k = 1:12
%             th_rad = deg2rad(all_theta_mapping(k)); % 获取当前位点的角度（弧度）
% 
%             % 判断当前位点属于阳极还是阴极
%             % 阳极
%             if ismember(k, anode_indices)
%                 % --- 使用阳极的最佳参数进行计算 ---
%                 % 检查是否满足穿孔条件：cos(θ) > cos(θc)
%                 if cos(th_rad) > cos(tc_a_rad)
%                     if j_a == 0
%                         % 当 j=0 时，公式简化为在穿孔区域内电导恒为 g0
%                         g_theta(k) = g0_a;
%                     else
%                         % 当 j>0 时，使用标准公式
%                         numerator = (abs(cos(th_rad)) - cos(tc_a_rad))^j_a;
%                         denominator = (1 - cos(tc_a_rad))^j_a;
%                         g_theta(k) = g0_a * numerator / denominator;
%                     end
%                 end
%             % 阴极
%             else % ismember(k, cathode_indices)
%                 % --- 使用阴极的最佳参数进行计算 ---
%                 % 检查是否满足穿孔条件：|cos(θ)| > cos(θc)
%                 if abs(cos(th_rad)) > cos(tc_c_rad)
%                     if j_c == 0
%                         % 当 j=0 时，公式简化为在穿孔区域内电导恒为 g0
%                         g_theta(k) = g0_c;
%                     else
%                         % 当 j>0 时，使用标准公式
%                         numerator = (abs(cos(th_rad)) - cos(tc_c_rad))^j_c;
%                         denominator = (1 - cos(tc_c_rad))^j_c;
%                         g_theta(k) = g0_c * numerator / denominator;
%                     end
%                 end
%             end
%         end
%         % 将计算出的一行12个点的电导率值存入最终结果矩阵
%         conductivity_data(t_idx, :) = g_theta;
%     end
% 
% 
%       %% 3. 绘图与结果可视化
%     % 此部分代码使用刚刚生成并保存在 conductivity_results_table 变量中的数据进行绘图
% 
%     fprintf('开始为 Sheet "%s" 生成结果图像...\n', sheet_name);
% 
%     % --- 提取绘图所需的数据 ---
%     time_vector = summary_table{:, 1}; % 更正：使用 summary_table 获取时间，它始终存在
%     conductivity_data_matrix = conductivity_data; % 更正：直接使用填充好的 conductivity_data 矩阵
%     num_time_points = length(time_vector);
%     num_angle_points = 12;
% 
%     % --- 图1: 电导率 vs. 角度位置 (时间作为不同曲线) ---
%     figure('Name', sprintf('Conductivity vs Angle for %s kV/cm', sheet_name), 'Visible', 'off'); % 创建一个不可见的figure窗口
%     hold on;
% 
%     colors = jet(num_time_points); % 使用jet颜色图为12条时间曲线生成不同颜色
%     legend_entries_time = cell(num_time_points, 1); % 初始化图例条目
% 
%     for t_idx = 1:num_time_points
%         % 绘制在 t_idx 时刻，12个角度位置的电导率
%         plot(1:num_angle_points, conductivity_data_matrix(t_idx, :), '-o', 'LineWidth', 1.5, 'Color', colors(t_idx,:));
%         legend_entries_time{t_idx} = sprintf('t = %.1f ns', time_vector(t_idx));
%     end
% 
%     hold off;
%     title(sprintf('Conductivity vs. Angle Location (%s kV/cm)', sheet_name));
%     xlabel('Angle Location Index (1-12)');
%     ylabel('Conductivity (S/m^2)'); % 您可以根据需要修改单位
%     xticks(1:num_angle_points); % 确保X轴刻度为整数
%     grid on;
%     legend(legend_entries_time, 'Location', 'northeastoutside'); % 将图例放在图的外面
% 
%     % 保存图1
%     plot1_filename = fullfile(output_folder, sprintf('Plot_Cond_vs_Angle_%skV.png', sheet_name));
%     saveas(gcf, plot1_filename); % gcf 获取当前figure句柄
%     fprintf('  - 成功保存图像: %s\n', plot1_filename);
% 
% 
%     % --- 图2: 电导率 vs. 时间 (角度作为不同曲线) ---
%     figure('Name', sprintf('Conductivity vs Time for %s kV/cm', sheet_name), 'Visible', 'off'); % 创建另一个不可见的figure窗口
%     hold on;
% 
%     colors = jet(num_angle_points); % 使用jet颜色图为12条角度曲线生成不同颜色
%     legend_entries_angle = cell(num_angle_points, 1);
% 
%     for angle_idx = 1:num_angle_points
%         % 【新增调试代码】打印出当前要绘制的角度和它的数据
%         fprintf('正在绘制角度 %d (Loc %d) 的数据:\n', all_theta_mapping(angle_idx), angle_idx);
%         disp(conductivity_data_matrix(:, angle_idx)'); % 使用'转置，让数据显示为一行
% 
%         % 绘制在 angle_idx 位置，所有时间点的电导率变化
%         plot(time_vector, conductivity_data_matrix(:, angle_idx), '-o', 'LineWidth', 1.5, 'Color', colors(angle_idx,:));
%         legend_entries_angle{angle_idx} = sprintf('Angle %d° (Loc %d)', all_theta_mapping(angle_idx), angle_idx);
%     end
% 
%     hold off;
%     title(sprintf('Conductivity vs. Time (%s kV/cm)', sheet_name));
%     xlabel('Time (ns)');
%     ylabel('Conductivity (S/m^2)'); % 您可以根据需要修改单位
%     grid on;
%     legend(legend_entries_angle, 'Location', 'northeastoutside');
% 
%     % 保存图2
%     plot2_filename = fullfile(output_folder, sprintf('Plot_Cond_vs_Time_%skV.png', sheet_name));
%     saveas(gcf, plot2_filename);
%     fprintf('  - 成功保存图像: %s\n', plot2_filename);
% 
%     close all; % 关闭所有打开的figure窗口，释放内存
% 
%     % --- 保存最终结果文件 ---
%     writetable(summary_table, param_summary_filename);
%     fprintf('成功写入参数摘要文件: %s\n', param_summary_filename);
%     conductivity_results_table = data(:, 1);
%     for k = 1:12, conductivity_results_table.(sprintf('Ori_%d_Conductivity', k)) = conductivity_data(:, k); end
%     writetable(conductivity_results_table, cond_filename);
%     fprintf('成功写入膜电导率结果文件: %s\n', cond_filename);
% 
% end
% fprintf('================================================\n');
% fprintf('所有电场条件的拟合和数据记录已全部完成。\n');