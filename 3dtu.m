
   clc; clear; close all;

    % ================= 1. 基础参数设置 =================
    L = 1;              % 梁长
    W = 0.2;            % 梁宽
    H = 0.02;           % 梁厚 (尽量薄，显得精致)
    N_points = 200;     % 【优化】加密网格，消除锯齿感
    max_amp = 0.4;      % 最大振幅 (端部位移量)
    
    % 【关键】定义5个位置：正最大，正中间，零，负中间，负最大
    % 顺序安排考虑绘图遮挡，通常先画中间再画两边，或者按顺序
    % 这里我们用一个数组定义所有状态
    amp_factors = [1.0, 0.5, 0, -0.5, -1.0]; 
    
    % 颜色风格 (冷灰色调，适合科研论文)
    beam_color = [0.75, 0.78, 0.82]; 
    fix_color = [0.85, 0.55, 0.4];   % 橙色固定端
    edge_color_solid = [0.2, 0.2, 0.2]; % 深灰色边框，比纯黑柔和

    % ================= 2. 窗口初始化 =================
    % 设置大尺寸窗口，提高预览清晰度
    figure('Color', 'w', 'Position', [100, 100, 1200, 800]);
    hold on; axis equal; axis off;
    
    % 设置视角 (等轴测偏侧视图)
    view(-40, 25); 

    % ================= 3. 计算几何数据 =================
    x = linspace(0, L, N_points);
    y = linspace(-W/2, W/2, 2); % 宽度方向只需要2个点
    [X, Y] = meshgrid(x, y);
    
    % 振型函数 (悬臂梁一阶近似: y = x^2 * (3L - x))
    % 归一化后乘以 L，保证末端位移为 L
    mode_shape = (x/L).^2 .* (3 - x/L) / 2; 

    % ================= 4. 循环绘制5个位置 =================
    for i = 1:length(amp_factors)
        factor = amp_factors(i);
        
        % 当前时刻的挠度 w(x)
        w = factor * max_amp * L * mode_shape;
        W_grid = repmat(w, 2, 1); % 扩展维度以匹配网格
        
        % --- 样式逻辑 (解决重叠显示差的核心) ---
        if abs(factor) == 1.0
            % 【极值位置 (+1, -1)】：最显眼
            face_alpha = 0.7;       % 较高不透明度
            current_edge = edge_color_solid; % 显示边框
            edge_alpha = 0.6;
        elseif factor == 0
            % 【零位置】：作为参考基准
            face_alpha = 0.15;      % 非常淡
            current_edge = 'none';  % 【关键】无边框，避免切割画面
            edge_alpha = 0;
        else
            % 【中间过渡位置 (+0.5, -0.5)】
            face_alpha = 0.3;       % 半透明
            current_edge = 'none';  % 【关键】无边框，消除重叠杂乱
            edge_alpha = 0;
        end
        
        % 绘制梁实体
        draw_beam_element(X, Y, W_grid, H, beam_color, ...
            face_alpha, current_edge, edge_alpha);
    end

    % ================= 5. 绘制固定端 =================
    draw_block([-0.1, 0], [-W/2, W/2], [-H*3, H*3], fix_color, 1.0);
    
    % ================= 6. 添加修饰 (力箭头与文字) =================
    % 在正负最大位置添加双向箭头示意
    quiver3(L*0.9, 0, max_amp*L*0.9 + 0.1, ...
            0, 0, -0.2, ... % 向下箭头
            0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(L*0.9, 0, -max_amp*L*0.9 - 0.1, ...
            0, 0, 0.2, ...  % 向上箭头
            0, 'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        
    text(L*1.05, 0, 0, 'Vibration Envelope', 'FontSize', 16, ...
        'FontName', 'Arial', 'HorizontalAlignment', 'left');

    % ================= 7. 光照与渲染优化 (解决质感差) =================
    lighting gouraud;    % 高洛德着色，使曲面光滑
    material dull;       % 哑光材质，减少刺眼反光
    
    % 多光源补光
    l1 = camlight('headlight');
    l2 = camlight('left');
    set(l1, 'Color', [1 1 1]*0.8); % 主光稍暗
    set(l2, 'Color', [1 1 1]*0.5); % 辅光柔和
    
    % 【关键】渲染器设置，保证透明度叠加正确
    set(gcf, 'Renderer', 'OpenGL'); 
    
    hold off;
    
    % ================= 8. 高清导出提示 =================
    fprintf('正在导出高清图片 (600 DPI)...\n');
    % 使用 exportgraphics 导出去白边的高清图
    exportgraphics(gcf, 'HighQuality_Vibration.png', 'Resolution', 600);
    fprintf('导出完成：HighQuality_Vibration.png\n');
end

% ================= 子函数：绘制梁单元 =================
function draw_beam_element(X, Y, W, h, color, f_alpha, e_color, e_alpha)
    % 计算上下表面
    Z_top = W + h/2;
    Z_bot = W - h/2;
    
    % 上表面
    surf(X, Y, Z_top, 'FaceColor', color, 'EdgeColor', e_color, ...
        'FaceAlpha', f_alpha, 'EdgeAlpha', e_alpha, 'MeshStyle', 'column');
    
    % 下表面
    surf(X, Y, Z_bot, 'FaceColor', color, 'EdgeColor', 'none', ...
        'FaceAlpha', f_alpha);
    
    % 侧面 (封闭几何体)
    % 提取边缘坐标构建侧面 patch
    x_edge = X(1,:); 
    w_edge = W(1,:); % 假设梁宽方向弯曲一致，取一行即可
    
    % 前侧面
    patch([x_edge, fliplr(x_edge)], ...
          [repmat(Y(1,1),1,length(x_edge)), repmat(Y(1,1),1,length(x_edge))], ...
          [w_edge+h/2, fliplr(w_edge-h/2)], ...
          color, 'EdgeColor', 'none', 'FaceAlpha', f_alpha);
      
    % 后侧面
    patch([x_edge, fliplr(x_edge)], ...
          [repmat(Y(2,1),1,length(x_edge)), repmat(Y(2,1),1,length(x_edge))], ...
          [w_edge+h/2, fliplr(w_edge-h/2)], ...
          color, 'EdgeColor', 'none', 'FaceAlpha', f_alpha);
      
    % 端面
    patch([X(1,end) X(2,end) X(2,end) X(1,end)], ...
          [Y(1,end) Y(2,end) Y(2,end) Y(1,end)], ...
          [Z_top(1,end) Z_top(2,end) Z_bot(2,end) Z_bot(1,end)], ...
          color, 'EdgeColor', e_color, 'FaceAlpha', f_alpha, 'EdgeAlpha', e_alpha);
end

% ================= 子函数：绘制固定块 =================
function draw_block(xr, yr, zr, color, alpha)
    V = [xr(1) yr(1) zr(1); xr(2) yr(1) zr(1); xr(2) yr(2) zr(1); xr(1) yr(2) zr(1);
         xr(1) yr(1) zr(2); xr(2) yr(1) zr(2); xr(2) yr(2) zr(2); xr(1) yr(2) zr(2)];
    F = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
    patch('Vertices', V, 'Faces', F, 'FaceColor', color, ...
        'FaceAlpha', alpha, 'EdgeColor', 'none');
end