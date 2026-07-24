function fitness = objective_thrust(x, k, m, ka, omega, v, disp_node, nodes_ext)
% 最大化推力的适应度函数
    V_front = x(1);
    V_rear = x(2);
    delta_phi = x(3);
    DC = x(4);
    
    % 构造4片压电片的电压向量
    omega0 = 2*pi*150;  % 激励频率
    t_eval = 0.5;  % 稳态时刻
    
    % 前段（P_UL, P_DL）
    V_UL = V_front * sin(omega0*t_eval);
    V_DL = -V_front * sin(omega0*t_eval) + DC;
    
    % 后段（P_UR, P_DR）
    V_UR = V_rear * sin(omega0*t_eval + delta_phi);
    V_DR = -V_rear * sin(omega0*t_eval + delta_phi) + DC;
    
    Va = [V_UL; V_UR; V_DL; V_DR];
    
    % 计算响应（频域快速求解）
    F = ka * Va;
    ac = 55.5; bc = 3.65e-06;
    C = ac*m + bc*k;
    A = k + 1i*omega0*C - omega0^2*m;
    q = A \ F;
    
    % 提取尾缘位移和速度
    [~, w] = ttd(disp_node, real(q), nodes_ext);
    idx_tail = find(nodes_ext(:,1) == max(nodes_ext(:,1)), 1);
    
    A_tail = abs(w(idx_tail));  % 尾缘振幅 (m)
    v_tail = omega0 * A_tail;   % 尾缘速度 (m/s)
    
    % 推力估计（简化 Lighthill 模型）
    thrust = 0.5 * 1000 * v_tail^2 * (0.01*0.001);  % ρ * v^2 * A
    
    % 惩罚项：避免 DC 过大
    penalty = 0;
    if abs(DC) > 250
        penalty = (abs(DC) - 250)^2 * 0.01;
    end
    
    % 适应度（GA最小化，所以返回负值）
    fitness = -thrust + penalty;
end