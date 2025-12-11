%% 计算跨膜电压和电导率
set(0,'defaultfigurecolor','w') ;
%% pre-processing 读取有效数据
% 该节用来对原始文件进行读取并保留有效数据
% 功能描述：
% 1.遍历每个xls原文件，使其包含每个sheet（时间点）成为一个N*12(N是该sheet下的样本量，12是位点数）的数据矩阵
% 2.计算该矩阵的均值和标准差分别放入最后2行,形成(N+2)*12的数据矩阵
% 3.生成对应的xls文件，保留在'output_file'路径下，后缀加上_processed

% 设置文件夹路径
folder_path = ''; % input path
output_folder = ''; % output path

% 确保输出文件夹存在
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
% 获取所有xls和xlsx文件
file_list = dir(fullfile(folder_path, '*.xls*'));
% 遍历所有文件
for i = 1:length(file_list)
    file_path = fullfile(folder_path, file_list(i).name);    
    % 获取所有表名
    sheets = sheetnames(file_path);    
    % 创建输出文件
    [~, file_name, ext] = fileparts(file_list(i).name);
    output_file = fullfile(output_folder, [file_name, '_processed', ext]);    
    % 遍历所有sheet
    for j = 1:length(sheets)
        % 读取数据
        data = readmatrix(file_path, 'Sheet', sheets(j));        
        % 数据处理
        if ~isempty(data)
            % 删除第二列包含NaN的行
            nan_indices = isnan(data(:,2));
            data(nan_indices,:) = [];          
            % 保留第三列及以后的数据
            if size(data, 2) >= 3
                data = data(:, 3:end);
            end      
            % 计算统计量
            if ~isempty(data)
                % 计算每列的均值
                mean_values = mean(data, 1);
                % 计算每列的标准差
                std_values = std(data, 0, 1); 
                % 添加均值和标准差行
                data = [data; mean_values; std_values];
                % 将处理后的数据写入新文件
                writematrix(data, output_file, 'Sheet', sheets(j));
            end
        end
    end
    
    fprintf('已处理文件: %s\n', file_list(i).name);
end

fprintf('所有文件预处理完成！\n');    
%% 批量计算各独立样本的跨膜电压及其统计数据
%
% 功能描述:
% 1. 遍历指定文件夹内的所有Excel文件 (*.xls, *.xlsx)。
% 2. 对每个文件，执行以下操作：
%    a. 读取文件中的每个sheet (时间点)。
%    b. 对每个独立样本数据计算其对应的跨膜电压值 (x1, x2)。
%    c. 以位点 (列) 为单位, 计算电压值的均值和标准差。
% 3. 将所有结果写入一个总的Excel文件中。每个输入文件会生成两个sheet：
%    - 一个以"_Mean"为后缀，包含电压均值。
%    - 一个以"_StdDev"为后缀，包含电压标准差。

%% 1. 初始化和参数设置
clear;         
clc;           
close all;     

% 输入文件夹路径 (包含所有result_*_processed.xls文件)
input_folder = ''; 

% 输出文件路径 (所有结果将保存在此文件中)
output_file = ''; 

% 跨膜电压-荧光强度高斯曲线参数 (根据 "映射结果" 提供)
a1 = 1.3438;
b1 = -552.8212;
c1 = 1.0062e+03;

%% 2. 批量数据处理
% 获取文件夹内所有xls和xlsx文件
file_list = dir(fullfile(input_folder, '*.xls*'));

% 检查是否找到了文件
if isempty(file_list)
    error('在指定文件夹中没有找到任何.xls或.xlsx文件: %s', input_folder);
end

% 删除旧的输出文件以防止写入错误和数据混淆
if isfile(output_file)
    delete(output_file);
    fprintf('已删除旧的输出文件: %s\n\n', output_file);
end

fprintf('共找到 %d 个文件需要处理。开始批量处理...\n', length(file_list));
fprintf('================================================\n');

% 生成通用的表头
num_sites = 12; % 位点数为12
headers_x1 = arrayfun(@(i) sprintf('位点%d_x1', i), 1:num_sites, 'UniformOutput', false);
headers_x2 = arrayfun(@(i) sprintf('位点%d_x2', i), 1:num_sites, 'UniformOutput', false);
headers_x1_std = arrayfun(@(i) sprintf('位点%d_x1_std', i), 1:num_sites, 'UniformOutput', false);
headers_x2_std = arrayfun(@(i) sprintf('位点%d_x2_std', i), 1:num_sites, 'UniformOutput', false);

header_mean = [{'Time(s)'}, headers_x1, headers_x2];
header_std = [{'Time(s)'}, headers_x1_std, headers_x2_std];

% 外层循环：遍历文件夹中的每一个文件
for i = 1:length(file_list)
    file_path = fullfile(input_folder, file_list(i).name);
    fprintf('--> 开始处理第 %d/%d 个文件: %s\n', i, length(file_list), file_list(i).name);
    
    % 获取源文件中的所有sheet名称 (每个sheet代表一个时间点)
    try
        sheets = sheetnames(file_path);
    catch ME
        fprintf('    [错误] 无法读取文件 "%s"，已跳过。原因: %s\n', file_list(i).name, ME.message);
        continue; % 跳过此文件，继续下一个
    end
    
    % 初始化用于存储当前文件结果的cell数组
    file_results_mean = header_mean;
    file_results_std = header_std;
    
    % 内层循环：遍历当前文件的每一个sheet (时间点)
    for s = 1:length(sheets)
        sheet_name = sheets{s};
        
        % 读取当前sheet的原始数据
        data = readmatrix(file_path, 'Sheet', sheet_name);
        
        if size(data, 1) <= 2
            fprintf('    - 警告: Sheet "%s" 数据行数不足，已跳过。\n', sheet_name);
            continue;
        end
        
        % 提取荧光强度数据 (排除最后两行的均值和标准差)
        fluo_data = data(1:end-2, :);
        
        if isempty(fluo_data)
            fprintf('    - 警告: Sheet "%s" 中没有有效的样本数据，已跳过。\n', sheet_name);
            continue;
        end
        
        % 通过单个数据计算每个sheet（时间点）的均值和标准差
        [n_samples, n_sites_actual] = size(fluo_data);
        voltages_x1 = NaN(n_samples, n_sites_actual);
        voltages_x2 = NaN(n_samples, n_sites_actual);
        
        for r = 1:n_samples
            for c = 1:n_sites_actual
                f = fluo_data(r, c);
                if ~isnan(f)
                    [x1, x2] = solve_voltage(f, a1, b1, c1); % 跨膜电压计算函数
                    voltages_x1(r, c) = abs(x1);
                    voltages_x2(r, c) = abs(x2);
                end
            end
        end
        
        % 计算均值和标准差
        mean_x1 = mean(voltages_x1, 1, 'omitnan');
        std_x1  = std(voltages_x1, 0, 1, 'omitnan');
        mean_x2 = mean(voltages_x2, 1, 'omitnan');
        std_x2  = std(voltages_x2, 0, 1, 'omitnan');
        
        % 将当前sheet的结果添加到总结果矩阵中
        current_row_mean = [mean_x1, mean_x2];
        current_row_std  = [std_x1, std_x2];
        
        file_results_mean(s+1, 1) = {sheet_name};
        file_results_mean(s+1, 2:end) = num2cell(current_row_mean);
        
        file_results_std(s+1, 1) = {sheet_name};
        file_results_std(s+1, 2:end) = num2cell(current_row_std);
    end
    
    % --- 将当前文件的结果写入输出文件 ---
    if size(file_results_mean, 1) > 1 % 确保有数据被处理
        % 生成sheet名称 (如 "result_5KV_processed" -> "5KV")
        [~, base_name, ~] = fileparts(file_list(i).name);
        sheet_base_name = regexprep(base_name, {'result_', '_processed'}, '');
        
        % 写入均值
        try
            writecell(file_results_mean, output_file, 'Sheet', [sheet_base_name, '_Mean']);
            fprintf('    + 成功写入均值数据到 Sheet: %s\n', [sheet_base_name, '_Mean']);
        catch ME
            fprintf('    [错误] 写入均值数据失败。原因: %s\n', ME.message);
        end
        
        % 写入标准差
        try
            writecell(file_results_std, output_file, 'Sheet', [sheet_base_name, '_StdDev']);
            fprintf('    + 成功写入标准差数据到 Sheet: %s\n', [sheet_base_name, '_StdDev']);
        catch ME
            fprintf('    [错误] 写入标准差数据失败。原因: %s\n', ME.message);
        end
    else
        fprintf('    - 文件 "%s" 未生成任何有效数据，已跳过写入步骤。\n', file_list(i).name);
    end
    fprintf('------------------------------------------------\n');
end

fprintf('\n批量处理完成！\n结果已保存至: %s\n', output_file);


%% 计算膜电导率 - temp
%
% 功能描述:
% 1. 自动读取输入Excel文件中的所有工作表。
% 2. 对每个工作表，使用两种方法计算12个位点的膜电导率。
%    - 方法1: 调用之前已定义的函数 calculate_sigma1。
%    - 方法2: 使用标准电穿孔公式。
% 3. 将每个工作表的计算结果保存到输出Excel文件中对应的同名工作表。
% 4. 输出格式为: Time, Ori 1(Method 1), Ori 1(Method 2), ...

clear; clc; close all;

% --- 1. 参数设置 ---

% 输入和输出文件
input_file = '';
output_file = '';

% 物理参数
R = 5e-6;             % 半径 (m)
d = 7e-9;             % 膜厚度 (m)
sigma_i = 0.5;          % 电导率 (S/m)
sigma_e = 1.0;          % 电导率 (S/m)

% 12个位点对应的角度 (单位: 度)
theta_degrees = [60, 30, 0, 30, 60, 90, 120, 150, 180, 150, 120, 90];

E_values_kV_per_cm = [5, 12, 50, 90, 100]; % 

% --- 2. 主程序 ---

fprintf('开始处理文件: %s\n', input_file);

% 获取输入Excel文件中的所有工作表名称
try
    sheets = sheetnames(input_file);
catch ME
    error('无法读取Excel文件，请检查文件名或文件路径是否正确。');
end

% 检查电场强度数组的数量是否与工作表数量匹配
if length(E_values_kV_per_cm) ~= length(sheets)
    warning('警告: 设置的电场强度数量(%d)与Excel工作表数量(%d)不匹配。程序将只处理前面部分。', ...
            length(E_values_kV_per_cm), length(sheets));
end

% 遍历每一个工作表
num_sheets_to_process = min(length(E_values_kV_per_cm), length(sheets));
for s = 1:num_sheets_to_process
    sheet_name = sheets{s};
    fprintf('\n正在处理工作表: %s\n', sheet_name);
    
    % --- 数据读取 ---
    data = readtable(input_file, 'Sheet', sheet_name);
    time_col = data.Time;
    n_sites = (width(data) - 1) / 2;
    
    % --- 计算准备 ---
    % 获取当前工作表对应的电场强度 E (V/m)
    current_E_kV_cm = E_values_kV_per_cm(s);
    E = current_E_kV_cm * 1e5; % 转换: 1 kV/cm = 1e5 V/m
    fprintf('  - 使用电场强度: %.1f kV/cm\n', current_E_kV_cm);

    % 初始化结果矩阵
    results_data = zeros(height(data), 2 * n_sites);
    
    % --- 遍历位点，调用两种方法计算 ---
    for site = 1:n_sites  
        % --- 方法2 (矢量化计算) ---
        % 一次性处理整列数据
        delta_V_vector = data{:, 2 * site} * 1e-3;
        delta_V0 = 1.5 * R * E * abs(cosd(theta_degrees(site)));
        numerator = 2 * sigma_e * sigma_i * d;
        denominator = R * (2 * sigma_e + sigma_i);
        sigma_m2_vector = (numerator / denominator) * ( 3*E*R./ 2*delta_V0 - 1);
        % ！！！
        results_data(:, 2 * site) = sigma_m2_vector;
        
        % --- 方法1 (逐点循环计算) ---
        % 这个方法必须对每个时间点的数据单独调用函数
        fprintf('    - 方法1: 正在逐点计算 (共 %d 个点), 请稍候...\n', height(data));
        for t = 1:height(data)
            delta_V_scalar = data{t, 2 * site} * 1e-3; % 获取单个数据点
            
            % 只有在数据点有效时才进行计算
            if ~isnan(delta_V_scalar) && delta_V_scalar > 0
                % 调用函数，传入单个数据点
                sigma_m1_scalar = calculate_sigma1(delta_V_scalar, E, theta_degrees(site), R, d, sigma_i, sigma_e);
                results_data(t, 2 * site - 1) = sigma_m1_scalar;
            end
        end
    end

    
    % --- 结果整理与输出 ---
    % 创建表头
    new_headers = cell(1, 2 * n_sites + 1);
    new_headers{1} = 'Time (ns)';
    for i = 1:n_sites
        new_headers{2*i}     = sprintf('Ori %d (Method 1)', i);
        new_headers{2*i + 1} = sprintf('Ori %d (Method 2)', i);
    end
    
    % 合并为表
    output_table = array2table([time_col, results_data]);
    output_table.Properties.VariableNames = new_headers;
    
    % 写入到输出Excel文件
    writetable(output_table, output_file, 'Sheet', sheet_name);
    fprintf('  - 工作表 "%s" 处理完成，结果已保存。\n', sheet_name);
end

fprintf('\n所有工作表处理完毕。\n');
fprintf('最终结果已保存至文件: %s\n', output_file);