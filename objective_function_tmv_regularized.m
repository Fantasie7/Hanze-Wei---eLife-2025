function [error, theory_tmv] = objective_function_tmv_regularized(params, exp_tmv_V, theta_deg, nu, N, j_param, pole_type, lambda)
%% 【正则化】拟合单个极的目标函数
    g0 = params(1);
    tc_rad = deg2rad(params(2));
    j = j_param;

    if strcmpi(pole_type, 'anode')
        g_func = @(th_rad) g0 * ((cos(th_rad) - cos(tc_rad)).^j) ./ ((1 - cos(tc_rad)).^j) .* (cos(th_rad) > cos(tc_rad));
    elseif strcmpi(pole_type, 'cathode')
        g_func = @(th_rad) g0 * ((abs(cos(th_rad)) - cos(tc_rad)).^j) ./ ((1 - cos(tc_rad)).^j) .* (abs(cos(th_rad)) > cos(tc_rad));
    else
        error('未知的pole_type');
    end

    odd_indices = 1:2:N;
    num_odd = length(odd_indices);
    g_mn = zeros(num_odd, num_odd);
    theta_int = linspace(0, pi, 200); x_int = cos(theta_int);
    for i = 1:num_odd
        m = odd_indices(i);
        for k = 1:num_odd
            n = odd_indices(k);
            g_theta_vals = g_func(theta_int);
            Pm = legendre_poly(m, x_int); 
            Pn = legendre_poly(n, x_int);
            integrand = Pm .* g_theta_vals .* Pn;
            g_mn(i, k) = trapz(x_int, integrand);
        end
    end

    A = zeros(num_odd, num_odd); 
    b = zeros(num_odd, 1);
    g_m1 = g_mn(:, 1);
    for i = 1:num_odd
        m = odd_indices(i);
        b(i) = (3/2) * g_m1(i);
        for k = 1:num_odd
            n = odd_indices(k);
            A(i, k) = (nu / (n + 1) + 1/n) * g_mn(i, k);
            if i == k
                A(i, k) = A(i, k) + 2 / (2*m + 1);
            end
        end
    end
    c_coeffs = A \ b;
    
    theory_tmv = zeros(size(theta_deg));
    x_exp = cosd(theta_deg);
    P1_val = legendre_poly(1, x_exp);
    theory_tmv = theory_tmv + (3/2) * P1_val;
    for i = 1:num_odd
        n = odd_indices(i);
        weight = -c_coeffs(i) * (1/n + nu/(n + 1));
        Pn_val = legendre_poly(n, x_exp);
        theory_tmv = theory_tmv + weight * Pn_val;
    end

    error_fit = sqrt(mean((abs(theory_tmv) - abs(exp_tmv_V)).^2));
    penalty_term = params(1)^2;
    error = error_fit + lambda * penalty_term;
    
    if ~isreal(error) || isnan(error), error = 1e9; end
end