clear; clc; close all;
%% 1. 参数定义 (单位: mm)
bm = 40.0;    % 致动器长度 (Actuator Length)
lm = 20.0;    % 致动器宽度 (Actuator Width)
be = 40.0;    % 鳍长度 (Fin Extension Length)
le = 55.0;    % 鳍尾端跨度 (Tail Span)
d_fork = 10.0; % 尾部分叉深度 (Fork Depth)

n_points = 100; % 采样密度

%% 2. 生成上边缘轮廓 (Upper Profile)
% --- A. 致动器段 (直线) ---
x_act = linspace(0, bm, n_points/2);
y_act = (lm / 2) * ones(size(x_act));

% --- B. 鳍面过渡段 (Smoothstep) ---
x_fin = linspace(bm, bm + be, n_points);
t = (x_fin - bm) / be; % 归一化坐标 [0, 1]

% 形状函数: Cubic Smoothstep (3t^2 - 2t^3)
% 在 t=0 处导数为0 (平滑连接直线)
% 在 t=1 处导数为0 (平滑到达最大宽度)
shape_func = 3*t.^2 - 2*t.^3; 

% 计算 Y 坐标
delta_width = (le - lm) / 2;
y_fin = (lm / 2) + delta_width * shape_func;

% --- C. 合并数据 ---
x_upper = [x_act, x_fin];
y_upper = [y_act, y_fin];

%% 3. 生成后缘轮廓 (Trailing Edge)
% 这是一个关于 y 的函数 x(y)
y_tail = linspace(le/2, -le/2, n_points);
y_norm = y_tail / (le/2); % 归一化 [-1, 1]

% 抛物线凹陷
x_tail = (bm + be) - d_fork * (1 - abs(y_norm).^2);

%% 4. 导出数据给 COMSOL
% COMSOL 的插值曲线通常接受两列数据 (x, y)
% 我们需要导出两段曲线：
% 1. 上边缘 (Upper Curve)
% 2. 后缘 (Trailing Edge)
data_upper = [x_upper', y_upper']; % 转置为列向量
data_tail  = [x_tail', y_tail'];   % 转置为列向量
% 保存为 txt 文件 (空格分隔)
writematrix(data_upper, 'Fin_Upper_Profile.txt', 'Delimiter', 'space');
writematrix(data_tail,  'Fin_Trailing_Edge.txt', 'Delimiter', 'space');

fprintf('数据已导出:\n 1. Fin_Upper_Profile.txt\n 2. Fin_Trailing_Edge.txt\n');

%% 5. 绘图验证
figure('Color', 'w');
hold on; axis equal; grid on;
% 绘制导出的曲线
plot(x_upper, y_upper, 'b-', 'LineWidth', 2);
plot(x_upper, -y_upper, 'b-', 'LineWidth', 2); % 镜像下边缘
plot(x_tail, y_tail, 'r-', 'LineWidth', 2);
% 封闭左侧
plot([0, 0], [-lm/2, lm/2], 'k-', 'LineWidth', 1);

title('COMSOL 几何轮廓预览');
xlabel('x (mm)'); ylabel('y (mm)');
legend('Upper Profile', 'Lower Profile', 'Trailing Edge');
hold off;





% =====================================================================
%  功能：读取 COMSOL .mphtxt 文件，仅提取四边形单元
%  适用：仅支持 4-节点壳单元的有限元程序
% =====================================================================
filename = '11.26wangge.mphtxt';
if ~isfile(filename)
    error('找不到文件: %s', filename);
end

fid = fopen(filename, 'r');
nodes = [];
q8_elems = []; % 存储8节点单元

fprintf('正在读取文件并转换单元类型...\n');

while ~feof(fid)
    line = fgetl(fid);

    % -----------------------------------------------------------------
    % 1. 读取节点坐标 (Vertices)
    % -----------------------------------------------------------------
    if contains(line, '# number of mesh vertices')
        num_nodes = sscanf(line, '%d');
        fgetl(fid); % 跳过 '0 # lowest mesh vertex index'
        fgetl(fid); % 跳过 空行
        fgetl(fid); % 跳过 '# Mesh vertex coordinates' 这一行文本
        % 读取坐标
        data = fscanf(fid, '%f %f %f', [3, num_nodes]);
        nodes = data';
        fprintf('  -> 已读取节点: %d 个\n', num_nodes);
    end

    % -----------------------------------------------------------------
    % 2. 读取 Quad2 (9节点) 并转换为 Quad8
    % -----------------------------------------------------------------
    % COMSOL mphtxt 中，二阶四边形标记为 "quad2"
    if contains(line, 'quad2 # type name')

        % 寻找元素数量
        while ~feof(fid)
            l = fgetl(fid);
            if contains(l, '# number of elements')
                num_elems = sscanf(l, '%d');
                break;
            end
        end

        % 跳过直到数据行
        while true
            pos = ftell(fid);
            l = fgetl(fid);
            if ~isempty(strtrim(l)) && ~startsWith(strtrim(l), '#')
                fseek(fid, pos, 'bof');
                break;
            end
        end

        % 读取 9 个节点 (COMSOL quad2 格式: 4角点 + 4边点 + 1中心点)
        % 我们读取 9 列，但在存储时只取前 8 列
        raw_data = fscanf(fid, '%d %d %d %d %d %d %d %d %d', [9, num_elems]);
        raw_data = raw_data';

        % 核心操作：丢弃第9列 (中心节点)，保留前8列
        q9_elems = raw_data(:, :);

        % 索引修正: 0-based -> 1-based
        q9_elems = q9_elems + 1;

        fprintf('  -> 已读取 9节点单元并转换为 8节点壳单元: %d 个\n', num_elems);
    end
end
fclose(fid);
if isempty(nodes) || isempty(q9_elems)
    error('读取失败：未找到节点或 quad2 单元数据。');
end
node_indices(1,1:8)= [1,2,4,3,5,8,9,6];
q8_elems(:,1:8) = q9_elems(:,node_indices);

% =====================================================================
%  3. 8节点单元可视化 (画出弯曲边缘并标注单元编号)
% =====================================================================
figure('Color', 'w', 'Name', '8-Node Shell Mesh with IDs');
hold on; axis equal; axis off;
title(['8-Node Serendipity Elements (Total: ', num2str(num_elems), ')']);

% 构建绘图顺序索引 (封闭环，用于绘制边界)
plot_order = [1, 5, 2, 6, 3, 7, 4, 8, 1];

% 循环绘制每个单元
for i = 1:num_elems
    % 1. 获取当前单元的节点索引
    idx_draw = q8_elems(i, plot_order); % 用于画线的顺序
    idx_all  = q8_elems(i, :);          % 该单元所有8个节点，用于计算形心
    
    % 2. 提取坐标
    xx = nodes(idx_draw, 1);
    yy = nodes(idx_draw, 2);
    
    % 3. 绘制单元边界
    plot(xx, yy, 'b-', 'LineWidth', 1.0);
    
    % 4. 计算单元几何中心 (作为标签位置)
    % 理论依据：x_c = (1/N) * sum(x_i)
    xc = mean(nodes(idx_all, 1));
    yc = mean(nodes(idx_all, 2));
    
    % 5. 标注单元编号
    % 'HorizontalAlignment', 'center' 保证文字居中于形心
    text(xc, yc, num2str(i), ...
        'Color', 'r', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end

% 调整视图
margin = 2;
xlim([min(nodes(:,1))-margin, max(nodes(:,1))+margin]);
ylim([min(nodes(:,2))-margin, max(nodes(:,2))+margin]);

fprintf('完成。已绘制单元网格并标注编号。\n');


% =====================================================================
%  3. 指定单元可视化 (画出弯曲边缘并标注单元编号)
% =====================================================================
% 设定要绘制的单元编号范围
draw_indices = 1:78; 

% 检查范围有效性
if max(draw_indices) > num_elems
    warning('请求的单元编号超出最大单元数，将自动截断。');
    draw_indices = draw_indices(draw_indices <= num_elems);
end

figure('Color', 'w', 'Name', 'Partial Mesh Visualization');
hold on; axis equal; axis off;
title(['Partial Mesh: Elements ' num2str(min(draw_indices)) ' to ' num2str(max(draw_indices))]);

% 构建绘图顺序索引 (封闭环)
plot_order = [1, 5, 2, 6, 3, 7, 4, 8, 1];

% 收集绘图涉及的节点索引，用于后续自动调整视图范围
active_node_indices = [];

% 循环绘制指定范围的单元
for i = draw_indices
    % 1. 获取数据
    idx_draw = q8_elems(i, plot_order); % 绘图连线顺序
    idx_all  = q8_elems(i, :);          % 计算形心用
    
    % 收集节点用于视图调整
    active_node_indices = [active_node_indices, idx_all]; %#ok<AGROW>
    
    xx = nodes(idx_draw, 1);
    yy = nodes(idx_draw, 2);
    
    % 2. 绘制边界
    plot(xx, yy, 'b-', 'LineWidth', 1.0);
    
    % 3. 计算形心并标注
    xc = mean(nodes(idx_all, 1));
    yc = mean(nodes(idx_all, 2));
    
    text(xc, yc, num2str(i), ...
        'Color', 'r', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end

% =====================================================================
%  视图自适应调整 (只聚焦于这78个单元)
% =====================================================================
active_node_indices = unique(active_node_indices);
if ~isempty(active_node_indices)
    active_nodes = nodes(active_node_indices, :);
    
    margin_ratio = 0.1; % 留白比例 10%
    x_range = range(active_nodes(:,1));
    y_range = range(active_nodes(:,2));
    
    xlim([min(active_nodes(:,1)) - x_range*margin_ratio, ...
          max(active_nodes(:,1)) + x_range*margin_ratio]);
      
    ylim([min(active_nodes(:,2)) - y_range*margin_ratio, ...
          max(active_nodes(:,2)) + y_range*margin_ratio]);
end

fprintf('完成。已绘制单元 %d 到 %d。\n', min(draw_indices), max(draw_indices));