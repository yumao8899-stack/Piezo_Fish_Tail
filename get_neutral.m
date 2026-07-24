function z0 = get_neutral(layers_data)
    % get_neutral_axis_offset: 
    % 计算非对称复合材料层合板的等效中性面相对于几何中心的偏移量 z₀。
    %
    % 输入:
    %   layers_data: 一个 N x 8 的矩阵，N是层数。每一行代表一层。
    %                行格式: [E1, E2, nu12, G12, G23, G13, thickness, angle_deg]
    %                - E1, E2: 主方向弹性模量 (e.g., in Pa)
    %                - nu12:   主泊松比
    %                - G12, G23, G13: 剪切模量 (e.g., in Pa)
    %                - thickness: 该层厚度 (e.g., in m)
    %                - angle_deg: 铺层角度 (单位: 度)
    % 输出:
    %   z0: 等效中性面相对于几何中心的偏移距离 (标量)。
    %       - z0 > 0: 中性面偏向 +z 方向 (上表面方向)。
    %       - z0 < 0: 中性面偏向 -z 方向 (下表面方向)。
    %       - z0 = 0: 中性面与几何中心重合 (对称铺层)。

    % 1. 计算总厚度 H 和每一层相对于几何中心的 z 坐标
    num_layers = size(layers_data, 1);
    thicknesses = layers_data(:, 7);
    total_thickness = sum(thicknesses);
    
    % z_coords(k) 是第 k-1 层的上表面/第 k 层的下表面坐标
    z_coords = zeros(num_layers + 1, 1);
    z_coords(1) = -total_thickness / 2; % 最底层的下表面坐标
    for k = 1:num_layers
        z_coords(k+1) = z_coords(k) + thicknesses(k);
    end

    % 2. 初始化 A (拉伸刚度) 和 B (拉伸-弯曲耦合) 矩阵
    A = zeros(3, 3);
    B = zeros(3, 3);

    % 3. 遍历每一层，计算其对 A 和 B 矩阵的贡献
    for k = 1:num_layers
        % 提取当前层的属性
        E1 = layers_data(k, 1);
        E2 = layers_data(k, 2);
        nu12 = layers_data(k, 3);
        G12 = layers_data(k, 4);
        angle_rad = layers_data(k, 8) * pi / 180;

        % 计算该层在材料主轴下的刚度矩阵 Q
        nu21 = nu12 * E2 / E1;
        Q11 = E1 / (1 - nu12 * nu21);
        Q22 = E2 / (1 - nu12 * nu21);
        Q12 = nu12 * E2 / (1 - nu12 * nu21);
        Q66 = G12;
        Q = [Q11, Q12, 0; Q12, Q22, 0; 0, 0, Q66];

        % 将 Q 矩阵旋转到全局坐标系，得到 Q_bar
        c = cos(angle_rad);
        s = sin(angle_rad);
        c2 = c^2;
        s2 = s^2;
        cs = c*s;
        
        Qxx = Q11*c2^2 + 2*(Q12+2*Q66)*s2*c2 + Q22*s2^2;
        Qyy = Q11*s2^2 + 2*(Q12+2*Q66)*s2*c2 + Q22*c2^2;
        Qxy = (Q11+Q22-4*Q66)*s2*c2 + Q12*(s2^2+c2^2);
        Qxs = (Q11-Q12-2*Q66)*s*c^3 + (Q12-Q22+2*Q66)*s^3*c;
        Qys = (Q11-Q12-2*Q66)*s^3*c + (Q12-Q22+2*Q66)*s*c^3;
        Qss = (Q11+Q22-2*Q12-2*Q66)*s2*c2 + Q66*(s2^2+c2^2);
        Q_bar = [Qxx, Qxy, Qxs; Qxy, Qyy, Qys; Qxs, Qys, Qss];
        
        % 获取当前层的上下表面坐标 (相对于几何中心)
        z_k_minus_1 = z_coords(k);
        z_k = z_coords(k+1);

        % 累加计算 A 和 B 矩阵
        A = A + Q_bar * (z_k - z_k_minus_1);
        B = B + (1/2) * Q_bar * (z_k^2 - z_k_minus_1^2);
    end
    
    % 4. 计算等效中性面偏移量 z0
    % B = A * z0_matrix  =>  z0_matrix = inv(A) * B
    % 我们使用一个等效的标量 z0 来近似
    if rcond(A) < 1e-15
        error('拉伸刚度矩阵 A 奇异或病态，无法计算中性面。请检查材料参数和铺层定义。');
    end
    z0_matrix = A \ B; % 使用左除 `\` 更稳定
    
    % 取对角线元素的平均值作为等效偏移量。
    % 对于大多数工程应用，这是一个合理的近似。
    z0 = mean(diag(z0_matrix)); 
    
end