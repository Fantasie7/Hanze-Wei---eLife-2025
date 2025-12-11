% 核心函数：用最小二乘求解膜电导率（允许非零残差）
function [sigma_m, residual] = calculate_sigma1(delta_V, E, theta_deg, R, d, sigma_i, sigma_e)
    theta = deg2rad(theta_deg);  % 角度转弧度
    
    % 定义目标函数：最小化 |eqn(sigma_m)|
    objective = @(sigma_m) abs( ...
        (3*sigma_e*(3*d*R^2*sigma_i + (3*d^2*R - d^3)*(sigma_m - sigma_i))) / ...
        (2*R^3*(sigma_m + 2*sigma_e)*(sigma_m + 0.5*sigma_i) - 2*((R-d)^3)*(sigma_e - sigma_m)*(sigma_i - sigma_m)) ...
        * E * R * abs(cos(theta)) - delta_V ...
    );
    
    % 最小二乘求解（设置搜索范围和优化选项）
    
    % 搜索范围：膜电导率通常为0到1000 S/m（可根据实际调整）
    lb = 1e-9;    % 下界（接近0，避免为负）
    ub = 1000;    % 上界
    
    % 初始猜测值（多组初始值提高稳定性）
    initial_guesses = [1e-6, 1e-3, 1, 10, 100];
    min_residual = Inf;
    best_sigma = NaN;
    
    % 尝试多组初始值，选择残差最小的解
    for guess = initial_guesses
        try
            [sigma, resnorm] = lsqnonlin(objective, guess, lb, ub);
            % 验证解的合理性（必须为正）
            if sigma > 0 && resnorm < min_residual
                min_residual = resnorm;
                best_sigma = sigma;
            end
        catch
            continue;  % 忽略求解失败的情况
        end
    end
    
    % 输出结果和残差
    sigma_m = best_sigma;
    residual = min_residual;
    
    % 过滤明显不合理的解（残差过大或电导率异常）
    if min_residual > 1e3 || (best_sigma <= 0 || best_sigma > 1000)
        sigma_m = NaN;
        residual = NaN;
    end
end