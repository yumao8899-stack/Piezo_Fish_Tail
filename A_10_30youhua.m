clc;clear
%% 1. 定义振型函数的参数
% --- x 方向的正弦波参数 ---
A1 = 0.8;      % 基础振幅 (最大值)
w1 = pi;       % x方向的角频率 (控制x方向有几个半波)
phi = 0;       % x方向的相位移
x0=1;

% --- y 方向的包络函数定义 ---
% 我们需要定义y的范围，因为包络通常与边界有关
y_min = -2;
y_max = 2;
L_y = y_max - y_min; % y方向的跨度

% 定义一个抛物线包络函数。
% 这个函数在 y = y_min 和 y = y_max 时为0，在 y 的中点处为 1。
% 这样的形状模拟了像吉他弦或简支梁的振动包络。
envelope = @(y) 1 - ((y - (y_min + L_y/2)) / (L_y/2)).^2;
% 你也可以尝试其他包络, 比如高斯包络:
% envelope = @(y) exp(-((y - (y_min + L_y/2)).^2) / 2);

% --- 组合成最终的振型函数 ---
target_shape1 = @(x,y) A1 * sin(w1 * x - phi) .* envelope(y);

vortex_control_shape = @(x,y) A1 * exp(-((x-x0)/pi).^2) .* sin(w1*y);
%% 2. 创建绘图网格
x_min = 0;
x_max = 2;
num_points_x = 100;
num_points_y = 100;
x_vec = linspace(x_min, x_max, num_points_x);
y_vec = linspace(y_min, y_max, num_points_y); 

% 使用 meshgrid 创建二维网格坐标
[X, Y] = meshgrid(x_vec, y_vec);
%% 3. 计算网格上每个点的 Z 值
Z = vortex_control_shape(X, Y);
%% 4. 绘制三维曲面图
figure('Name', 'Vibration Mode Shape with Envelope', 'NumberTitle', 'off');
surf(X, Y, Z);
% --- 图形美化 ---
title('Vibration Shape with Y-direction Envelope');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Amplitude (Z)');
colorbar;
shading interp;
view(35, 25); % 调整视角
grid on;
axis tight; % 让图形充满坐标轴