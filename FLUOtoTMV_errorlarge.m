clc;
clear all;
%% 输入数据
X = [-3600, -3400, -3200, -3000, -2800, -2600, -2400, -2200, -2000, -1800,...
     -1600, -1400, -1200, -1000, -800, -600, -400, -300, -200, -100,...
     0, 100, 200, 300, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]';
y = [2.01e-5, 9.06e-5, 4.43e-4, 1.78e-3, 5.84e-3, 1.64e-2, 4.03e-2,...
     8.68e-2, 1.67e-1, 2.91e-1, 4.62e-1, 6.69e-1, 8.94e-1, 1.10, 1.19,...
     1.26, 1.31, 1.33, 1.26, 1.10, 0.89, 0.67, 0.46, 0.29, 0.17, 0.087,...
     0.040, 0.016, 0.0058, 0.0018, 0.00054, 0.00013, 2.72e-4]';

%% 设置拟合选项和初始参数
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
                      'Display', 'off', 'MaxIterations', 2000, ...
                      'FunctionTolerance', 1e-9, 'StepTolerance', 1e-9);

% 三个高斯函数叠加模型 (9个参数)
model = @(p,x) p(1)*exp(-((x-p(2))/p(3)).^2) + ...
               p(4)*exp(-((x-p(5))/p(6)).^2) + ...
               p(7)*exp(-((x-p(8))/p(9)).^2);

% 初始参数估计 (基于数据特征)
p0 = [1.3, -300, 100,    % 主峰 (中心在-300附近)
      0.3, -1000, 500,   % 左侧肩部
      0.2, 200, 300];     % 右侧尾部

% 参数边界 (振幅和宽度为正)
lb = [0, -inf, 0, 0, -inf, 0, 0, -inf, 0];
ub = [];

%% 加权最小二乘拟合 (权重=1/y)
weight = 1./y;
objective = @(p) weight .* (model(p, X) - y);

% 非线性拟合
[p_opt, ~, residual, ~, ~] = lsqnonlin(objective, p0, lb, ub, options);

%% 计算拟合值和相对误差
y_fit = model(p_opt, X);
rel_error = abs(y_fit - y) ./ y;
max_rel_error = max(rel_error);

%% 验证拟合精度
fprintf('最大相对误差: %.4f%%\n', max_rel_error*100);
if max_rel_error < 0.001
    fprintf('满足0.1%%精度要求\n');
else
    fprintf('警告: 未达到0.1%%精度要求\n');
end

%% 可视化结果
figure('Position', [100, 100, 1200, 500])
subplot(1,2,1)
semilogy(X, y, 'bo', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', '原始数据')
hold on
semilogy(X, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', '拟合曲线')
xlabel('X')
ylabel('y (对数刻度)')
title('原始数据与拟合曲线')
legend('Location', 'best')
grid on

subplot(1,2,2)
plot(X, rel_error*100, 'ks-', 'LineWidth', 1.5, 'MarkerFaceColor', 'g')
xlabel('X')
ylabel('相对误差 (%)')
title('拟合相对误差')
yline(0.1, 'r--', '0.1% 阈值', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left')
grid on
ylim([0, max(0.2, 1.1*max(rel_error*100))])

%% 输出拟合函数
fprintf('\n拟合函数表达式:\n');
fprintf('y = %.6f * exp(-((x - (%.6f))/%.6f)^2) + ...\n', p_opt(1:3));
fprintf('     %.6f * exp(-((x - (%.6f))/%.6f)^2) + ...\n', p_opt(4:6));
fprintf('     %.6f * exp(-((x - (%.6f))/%.6f)^2)\n', p_opt(7:9));

%% 高精度评估函数
fit_func = @(x) p_opt(1)*exp(-((x-p_opt(2))/p_opt(3)).^2) + ...
                p_opt(4)*exp(-((x-p_opt(5))/p_opt(6)).^2) + ...
                p_opt(7)*exp(-((x-p_opt(8))/p_opt(9)).^2);

% 测试拟合函数
test_x = -3500:100:2500;
test_y = fit_func(test_x);

% 绘制完整拟合曲线
figure
semilogy(X, y, 'bo', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', '原始数据')
hold on
semilogy(test_x, test_y, 'r-', 'LineWidth', 1.5, 'DisplayName', '拟合函数')
xlabel('X')
ylabel('y (对数刻度)')
title('高精度拟合函数')
legend('Location', 'best')
grid on