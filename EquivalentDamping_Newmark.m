function [t_out, X_out] = EquivalentDamping_Newmark(M, K_complex, F_matrix, t, f_ex)
% SOLVE_EQUIVALENTDAMPING_NEWMARK 使用复模量等效阻尼和 Newmark-beta 法求解动力学响应
% 
% 输入参数:
%   M         - 整体质量矩阵 (N x N)
%   K_complex - 整体复刚度矩阵 (N x N), 由复模量组装而来
%   F_matrix  - 随时间变化的外力矩阵 (N x length(t)), 每一列对应一个时间步的节点力
%   t         - 时间向量 (1 x N_steps)
%   f_ex      - 结构当前的主导激励频率 (Hz), 用于换算等效黏性阻尼
%
% 输出参数:
%   t_out     - 时间向量
%   X_out     - 节点位移矩阵 (N x length(t))

    fprintf('====== 开始基于等效黏性阻尼的时域积分 ======\n');
    
    %% 1. 核心阻尼转换 (严格对应论文逻辑)
    % 提取主导圆频率
    omega_0 = 2 * pi * f_ex; 
    
    % 提取实部刚度与虚部阻尼
    K = real(K_complex);
    C = imag(K_complex) / omega_0; 
    
    % 检查矩阵是否对称（有限元矩阵通常应保持对称性）
    C = (C + C') / 2; 
    K = (K + K') / 2;
    
    %% 2. Newmark-beta 算法初始化 (无条件稳定格式)
    % 参数设定
    gamma = 0.5;
    beta  = 0.25;
    
    dt = t(2) - t(1);
    N_dof = size(M, 1);
    N_steps = length(t);
    
    % 初始化位移、速度、加速度矩阵
    X = zeros(N_dof, N_steps);
    V = zeros(N_dof, N_steps);
    A = zeros(N_dof, N_steps);
    
    % 初始状态计算 (假设初始位移和速度为 0)
    % M*A_0 = F_0 - C*V_0 - K*X_0
    A(:, 1) = M \ (F_matrix(:, 1)); 
    
    %% 3. 计算等效刚度矩阵并进行分解 (提升循环内求解效率)
    % K_hat = K + (gamma/(beta*dt))*C + (1/(beta*dt^2))*M
    a0 = 1 / (beta * dt^2);
    a1 = gamma / (beta * dt);
    a2 = 1 / (beta * dt);
    a3 = 1 / (2 * beta) - 1;
    a4 = gamma / beta - 1;
    a5 = (dt / 2) * (gamma / beta - 2);
    a6 = dt * (1 - gamma);
    a7 = gamma * dt;
    
    K_hat = K + a1 * C + a0 * M;
    
    % 提前对等效刚度矩阵进行 Cholesky 或 LU 分解以加速运算
    fprintf('正在分解等效刚度矩阵...\n');
    [L, U, P] = lu(K_hat);
    
    %% 4. 逐步时间积分
    fprintf('正在进行时域 Newmark-beta 积分, 总步数: %d...\n', N_steps);
    
    for i = 1:(N_steps - 1)
        % 计算等效载荷向量 P_hat
        % P_hat = F_(t+dt) + M*[a0*X_t + a2*V_t + a3*A_t] + C*[a1*X_t + a4*V_t + a5*A_t]
        
        vec_M = a0 * X(:, i) + a2 * V(:, i) + a3 * A(:, i);
        vec_C = a1 * X(:, i) + a4 * V(:, i) + a5 * A(:, i);
        
        P_hat = F_matrix(:, i+1) + M * vec_M + C * vec_C;
        
        % 解代数方程求 t+dt 时刻的位移: K_hat * X_(t+dt) = P_hat
        % 使用预先分解的 L U P 极速求解
        X(:, i+1) = U \ (L \ (P * P_hat));
        
        % 更新 t+dt 时刻的加速度和速度
        A(:, i+1) = a0 * (X(:, i+1) - X(:, i)) - a2 * V(:, i) - a3 * A(:, i);
        V(:, i+1) = V(:, i) + a6 * A(:, i) + a7 * A(:, i+1);
        
        % 进度条提示 (每 10% 打印一次)
        if mod(i, floor(N_steps/10)) == 0
            fprintf('  完成进度: %d%%\n', round((i/N_steps)*100));
        end
    end
    
    t_out = t;
    X_out = X;
    fprintf('====== 时域积分完成 ======\n');
end