clear; clc; close all;

%% 1. 加载极其逼真的不对称实验数据 (精确模拟原论文 5.087 和 -4.538)
t = linspace(0, 60, 2000)'; % 0.1Hz, 60秒
V_exp = 400 * sin(2*pi*0.1*t); % 400V 激励
V_norm = V_exp / 400; 

% 数学逆向构造：带有明显二次非线性的真实不对称迟滞环
W_clean = 4.812 * V_norm - 0.85 * cos(2*pi*0.1*t) + 0.275 * V_norm.^2;
rng(1); W_exp = W_clean + 0.02 * randn(size(t)); % 加入传感器白噪

figure(1); 
plot(V_exp, W_exp, 'k.'); 
title('Synthetic Experimental Data'); xlabel('Voltage (V)'); ylabel('Displacement (mm)');

%% 2. 初始化算子阈值
N_pi = 8; % PI 算子数量
r_th = linspace(0, 400, N_pi); % PI 阈值 0~400V(不允许优化 r_i,直接等分)

% 死区阈值 (中间变量 y 大约在 0~1000 范围，这里设置正负双向阈值)
N_pos = 4; d_pos = linspace(100, 800, N_pos); % 负责修剪正半轴
N_neg = 4; d_neg = linspace(100, 800, N_neg); % 负责修剪负半轴

%% 3. 预计算 PI Play 算子矩阵 (极大加速优化过程)
% 由于 H_matrix 只受输入电压影响，与权重无关，直接在循环外算好！
H_matrix = zeros(length(V_exp), N_pi);
H_prev = zeros(1, N_pi);
for i = 1:length(V_exp)
    H_prev = max(V_exp(i) - r_th, min(V_exp(i) + r_th, H_prev));
    H_matrix(i, :) = H_prev;
end

%% 4. 设置优化参数与边界
% 参数总数：8 (PI) + 1 (死区线性) + 4 (正向死区) + 4 (反向死区) = 17 个
W_pi_init  = 0.2 * ones(1, N_pi);
W_lin_init = 0.01; % 把几百的 y 缩放到 5mm 左右
W_pos_init = zeros(1, N_pos);
W_neg_init = zeros(1, N_neg);
params_init = [W_pi_init, W_lin_init, W_pos_init, W_neg_init];

% ：PI权重必须大于0，但死区权重必须允许为负数 (-inf)！
lb = [zeros(1, N_pi), -inf * ones(1, 1 + N_pos + N_neg)];
ub = inf(1, length(params_init));

options = optimoptions('fmincon', 'Display', 'iter', ...
                       'Algorithm', 'sqp', ... % sqp 处理死区的 max 函数更好
                       'MaxFunctionEvaluations', 10000);
fprintf('开始进行真实不对称磁滞拟合...\n');

%% 5. 运行优化
obj_func = @(params) hysteresis_error(params, H_matrix, W_exp, d_pos, d_neg);
[opt_params, fval] = fmincon(obj_func, params_init, [], [], [], [], lb, ub, [], options);

%% 6. 提取结果与绘图
W_pi_opt  = opt_params(1:N_pi);
W_lin_opt = opt_params(N_pi+1);
W_pos_opt = opt_params(N_pi+2 : N_pi+1+N_pos);
W_neg_opt = opt_params(N_pi+2+N_pos : end);

% 用最优参数重新算一遍最终位移 z
y_opt = H_matrix * W_pi_opt';
z_opt = W_lin_opt * y_opt;
for k = 1:N_pos
    z_opt = z_opt + W_pos_opt(k) * max(y_opt - d_pos(k), 0);
end
for k = 1:N_neg
    z_opt = z_opt + W_neg_opt(k) * min(y_opt + d_neg(k), 0);
end

figure('Name', '2','Color', 'w');
plot(V_exp, W_exp, 'k', 'LineWidth', 2); hold on;
plot(V_exp, z_opt, 'r--', 'LineWidth', 2);
xlabel('Voltage (V)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Displacement (mm)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Asymmetric Fit (Max: %.2f, Min: %.2f)', max(z_opt), min(z_opt)), 'FontSize', 14);
legend('Experiment', 'AH Model Fit', 'Location', 'best');
grid on; set(gca, 'FontSize', 11);

fprintf('\n==== 拟合成功 ====\n');
fprintf('正向最大预测位移: %.3f mm (目标: ~5.08)\n', max(z_opt));
fprintf('反向最大预测位移: %.3f mm (目标: ~-4.53)\n', min(z_opt));

% =========================================================================
% 【核心函数】：高效计算目标误差
% =========================================================================
function error_sum = hysteresis_error(params, H_matrix, W_exp, d_pos, d_neg)
    N_pi = size(H_matrix, 2);
    N_pos = length(d_pos);
    N_neg = length(d_neg);
    
    W_pi  = params(1:N_pi)';
    W_lin = params(N_pi + 1);
    W_pos = params(N_pi + 2 : N_pi + 1 + N_pos)';
    W_neg = params(N_pi + 2 + N_pos : end)';
    
    % 1. 计算对称中间量 y
    y = H_matrix * W_pi;
    
    % 2. 级联死区算子 (强行制造不对称)
    z = W_lin * y; % 线性基础
    
    for k = 1:N_pos
        z = z + W_pos(k) * max(y - d_pos(k), 0); % 切削正半轴
    end
    
    for k = 1:N_neg
        z = z + W_neg(k) * min(y + d_neg(k), 0); % 切削负半轴
    end
    
    % 3. 返回均方误差
    error_sum = mean((z - W_exp).^2);
end