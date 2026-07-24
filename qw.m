% =========================================================
% 【诊断工具】：手动扫描 a1，看看振幅到底变不变！
% =========================================================
sweep_points = 40;
a1_range = linspace(-30, 30, sweep_points);
a2_range = linspace(-30, 30, sweep_points);
b1_range = linspace(-10, 10, sweep_points); % 阻尼通常较难出现极端的负数

% 初始化存储数组
amp_a1 = zeros(1, sweep_points);
amp_a2 = zeros(1, sweep_points);
amp_b1 = zeros(1, sweep_points);
disp('正在进行手动参数扫描诊断，请稍候...');
 omega_test = 2 * pi * 7;


fprintf('正在扫描 a1...\n');
for idx = 1:sweep_points
    Mw_temp = zeros(sys_dof, sys_dof); Cw_temp = zeros(sys_dof, sys_dof);
    for loopi = 1:size(elements, 1)
        zmtemp_local = zmtemp_cell{loopi};

        % 调用水动力基底矩阵
        [mw_eff, cw_eff] = wateradding(zmtemp_local, rho_fluid, mu, omega_test, a1_range(idx), 1, 1);
        % 节点映射组装
        index = zeros(40, 1);
        for zi=1:8
            index((zi-1)*5 + 1) = disp_node(elements(loopi, zi), 1);
            index((zi-1)*5 + 2) = disp_node(elements(loopi, zi), 2);
            index((zi-1)*5 + 3) = disp_node(elements(loopi, zi), 3);
            index((zi-1)*5 + 4) = disp_node(elements(loopi, zi), 4);
            index((zi-1)*5 + 5) = disp_node(elements(loopi, zi), 5);
        end

        nonzero_idx = index(index ~= 0);
        valid_mask = index ~= 0;

        Mw_global(nonzero_idx, nonzero_idx) = Mw_global(nonzero_idx, nonzero_idx) + mw_eff(valid_mask, valid_mask);
        Cw_global(nonzero_idx, nonzero_idx) = Cw_global(nonzero_idx, nonzero_idx) + cw_eff(valid_mask, valid_mask);
    end
    M_total = m + Mw_global;
    C_total = C_matrix + Cw_global; % 结构瑞利阻尼与流体阻尼直接叠加

    % 标准频域代数方程: (-w^2 M + iw C + K) * U = F
    H_dyn = k - (omega_test^2) * M_total + (1i * omega_test) * C_total;
    U = H_dyn \ F_amp_vector;
    amp_a1(idx) = abs(U(jth)) * 1000;
end

fprintf('正在扫描 a2...\n');
for idx = 1:sweep_points
    Mw_temp = zeros(sys_dof, sys_dof); Cw_temp = zeros(sys_dof, sys_dof);
    for loopi = 1:size(elements, 1)
        [mw_eff, cw_eff] = wateradding(zmtemp_local, rho_fluid, mu, omega_test, 1, a2_range(idx), 1);
        % 节点映射组装
        index = zeros(40, 1);
        for zi=1:8
            index((zi-1)*5 + 1) = disp_node(elements(loopi, zi), 1);
            index((zi-1)*5 + 2) = disp_node(elements(loopi, zi), 2);
            index((zi-1)*5 + 3) = disp_node(elements(loopi, zi), 3);
            index((zi-1)*5 + 4) = disp_node(elements(loopi, zi), 4);
            index((zi-1)*5 + 5) = disp_node(elements(loopi, zi), 5);
        end

        nonzero_idx = index(index ~= 0);
        valid_mask = index ~= 0;

        Mw_global(nonzero_idx, nonzero_idx) = Mw_global(nonzero_idx, nonzero_idx) + mw_eff(valid_mask, valid_mask);
        Cw_global(nonzero_idx, nonzero_idx) = Cw_global(nonzero_idx, nonzero_idx) + cw_eff(valid_mask, valid_mask);
    end
    M_total = m + Mw_global;
    C_total = C_matrix + Cw_global; % 结构瑞利阻尼与流体阻尼直接叠加
    
    H_dyn = k - (omega_test^2)*(m + Mw_temp) + 1i*omega_test*(C_matrix + Cw_temp);
    U_complex = H_dyn \ F_amp_vector;
    amp_a2(idx) = abs(U(jth)) * 1000;
end

% ---------------------------------------------------------
% 3. 扫描 b1 (固定 a1=0, a2=0)
% ---------------------------------------------------------
fprintf('正在扫描 b1...\n');
for idx = 1:sweep_points
    Mw_temp = zeros(sys_dof, sys_dof); Cw_temp = zeros(sys_dof, sys_dof);
    for loopi = 1:size(elements, 1)
        [mw_eff, cw_eff] = wateradding(zmtemp_local, rho_fluid, mu, omega_test, 1, 1, b1_range(idx));
        

        % 节点映射组装
        index = zeros(40, 1);
        for zi=1:8
            index((zi-1)*5 + 1) = disp_node(elements(loopi, zi), 1);
            index((zi-1)*5 + 2) = disp_node(elements(loopi, zi), 2);
            index((zi-1)*5 + 3) = disp_node(elements(loopi, zi), 3);
            index((zi-1)*5 + 4) = disp_node(elements(loopi, zi), 4);
            index((zi-1)*5 + 5) = disp_node(elements(loopi, zi), 5);
        end

        nonzero_idx = index(index ~= 0);
        valid_mask = index ~= 0;

        Mw_global(nonzero_idx, nonzero_idx) = Mw_global(nonzero_idx, nonzero_idx) + mw_eff(valid_mask, valid_mask);
        Cw_global(nonzero_idx, nonzero_idx) = Cw_global(nonzero_idx, nonzero_idx) + cw_eff(valid_mask, valid_mask);
    end
    M_total = m + Mw_global;
    C_total = C_matrix + Cw_global; % 结构瑞利阻尼与流体阻尼直接叠加
   
    H_dyn = k - (omega_test^2)*(m + Mw_temp) + 1i*omega_test*(C_matrix + Cw_temp);
    U_complex = H_dyn \ F_amp_vector;
    amp_b1(idx) = abs(U_complex(jth)) * 1000;
end

% =========================================================================
% 绘制 1x3 敏感度对比图
% =========================================================================
figure('Name', '参数敏感度扫描诊断', 'Position', [100, 100, 1200, 400]);

subplot(1,3,1);
plot(a1_range, amp_a1, 'b-o', 'LineWidth', 1.5);
xlabel('a_1 (附加质量修正)'); ylabel('稳态振幅 (mm)');
title('a_1 敏感度曲线'); grid on;
yline(3.0, 'r--', '实验目标 3mm', 'LineWidth', 1.5); % 标出你的实验靶点

subplot(1,3,2);
plot(a2_range, amp_a2, 'g-o', 'LineWidth', 1.5);
xlabel('a_2 (附加质量修正)');
title('a_2 敏感度曲线'); grid on;
yline(3.0, 'r--', 'LineWidth', 1.5);

subplot(1,3,3);
plot(b1_range, amp_b1, 'k-o', 'LineWidth', 1.5);
xlabel('b_1 (附加阻尼修正)');
title('b_1 敏感度曲线'); grid on;
yline(3.0, 'r--', 'LineWidth', 1.5);


