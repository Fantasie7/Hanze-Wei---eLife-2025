function [error, theory_tmv_V] = objective_function_tmv(p, exp_tmv_V, theta_deg, ...
    E0_V_per_cm, a_cm, kappa_i_S_per_cm, nu, N, j_val, pole_type)
%% 目标函数：基于H.1991无量纲化框架
% 输入：
%   p: [实际G0 (S/cm²), theta_c (deg)]
%   exp_tmv_V: 实验实际跨膜电位（V）
%   E0_V_per_cm: 当前实验E0（V/cm）
%   a_cm: 细胞半径（cm），kappa_i_S_per_cm: 细胞内电导率（S/cm）
% 输出：
%   error: 实际跨膜电位的均方根误差
%   theory_tmv_V: 理论实际跨膜电位（V）

    % 1. 解包参数（实际值）
    G0 = p(1);                % 实际膜电导（S/cm²）
    tc_deg = p(2);            % 临界角度（deg）
    tc_rad = deg2rad(tc_deg);
    j = j_val;

    % 2. 无量纲化转换（实际值→文献无量纲变量）
    % 2.1 无量纲膜电导g0（文献Eq.A8）
    g0 = (a_cm / kappa_i_S_per_cm) * G0;  
    % 2.2 实验无量纲跨膜电位ψ_exp（文献Eq.A5）
    psi_exp = exp_tmv_V / (a_cm * E0_V_per_cm);  

    % 3. 文献膜电导分布g(θ)（无量纲，Eq.A16）
    if strcmpi(pole_type, 'anode')
        g_func = @(th_rad) g0 * ((abs(cos(th_rad)) - cos(tc_rad)).^j) ./ ...
            ((1 - cos(tc_rad)).^j) .* (abs(cos(th_rad)) > cos(tc_rad));
    elseif strcmpi(pole_type, 'cathode')
        g_func = @(th_rad) g0 * ((abs(cos(th_rad)) - cos(tc_rad)).^j) ./ ...
            ((1 - cos(tc_rad)).^j) .* (abs(cos(th_rad)) > cos(tc_rad));
    else
        error('pole_type必须为"anode"或"cathode"（Hibino et al. 1991）');
    end

    % 4. 文献线性方程组求解（无量纲ψ，Eq.A13）
    % 4.1 计算g_mn矩阵（无量纲）
    odd_indices = 1:2:N;
    num_odd = length(odd_indices);
    g_mn = zeros(num_odd, num_odd);
    theta_int = linspace(0, pi, 200);  % 加密积分点，提高精度
    x_int = cos(theta_int);
    for i = 1:num_odd
        m = odd_indices(i);
        for k = 1:num_odd
            n = odd_indices(k);
            g_theta_vals = g_func(theta_int);
            Pm = legendre_poly(m, x_int);
            Pn = legendre_poly(n, x_int);
            integrand = Pm .* g_theta_vals .* Pn;
            g_mn(i, k) = trapz(x_int, integrand);  % 数值积分（文献附录方法）
        end
    end

    % 4.2 构建A*c = b（无量纲，Eq.A13）
    A = zeros(num_odd, num_odd);
    b = zeros(num_odd, 1);
    g_m1 = g_mn(:, 1);
    for i = 1:num_odd
        m = odd_indices(i);
        b(i) = (3/2) * g_m1(i);  % 文献b向量定义
        for k = 1:num_odd
            n = odd_indices(k);
            A(i, k) = (nu/(n+1) + 1/n) * g_mn(i, k);
            if i == k
                A(i, k) = A(i, k) + 2/(2*m + 1);  % 克罗内克函数项
            end
        end
    end
    c_coeffs = A \ b;  % 求解系数c（无量纲）

    % 5. 计算理论无量纲跨膜电位ψ_theory（文献Eq.A10推导）
    x_exp = cosd(theta_deg);
    psi_theory = (3/2) * legendre_poly(1, x_exp);  % 基础项（n=1）
    for i = 2:num_odd
        n = odd_indices(i);
        weight = -c_coeffs(i) * (1/n + nu/(n+1));
        psi_theory = psi_theory + weight * legendre_poly(n, x_exp);  % 求和项
    end

    % 6. 无量纲→实际值转换（ψ→ΔΨ，文献Eq.A5逆运算）
    theory_tmv_V = psi_theory * (a_cm * E0_V_per_cm);  % 理论实际跨膜电位（V）

    % 7. 误差计算（实际值对比，保留符号，符合文献荧光-电位线性关系，🔶1-46）
    error = sqrt(mean((theory_tmv_V - exp_tmv_V).^2));
    % 异常值处理
    if ~isreal(error) || isnan(error) || error > 1e3
        error = 1e9;
    end
end