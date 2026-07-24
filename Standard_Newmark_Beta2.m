function [t, X] = Standard_Newmark_Beta2(M, C, K, F, t)
    % Standard_Newmark_Beta: 实数域 Newmark-beta 时域积分 (修正版)
    % 输入:
    %   M, C, K - 系统的质量、阻尼、实数刚度矩阵
    %   F       - 激振力矩阵 (每一列对应一个时间步的力向量)
    %   t       - 时间向量
    % 输出:
    %   t       - 时间向量
    %   X       - 位移响应矩阵 (每一列对应一个时间步的位移向量)

    % Newmark-beta 算法参数 (平均加速度法，无条件稳定)
    gamma = 0.5;
    beta  = 0.25;

    dt = t(2) - t(1);
    Nt = length(t);
    N_dof = size(M, 1);

    % 初始化响应矩阵
    X = zeros(N_dof, Nt);    % 位移
    V = zeros(N_dof, Nt);    % 速度
    A = zeros(N_dof, Nt);    % 加速度

    % 初始条件 (静止起振)
    X(:,1) = zeros(N_dof, 1);
    V(:,1) = zeros(N_dof, 1);
    A(:,1) = M \ (F(:,1) - C*V(:,1) - K*X(:,1));

    % 积分常数预计算
    a0 = 1 / (beta * dt^2);
    a1 = gamma / (beta * dt);
    a2 = 1 / (beta * dt);
    a3 = 1 / (2*beta) - 1;
    a4 = gamma / beta - 1;
    a5 = (dt / 2) * (gamma / beta - 2);

    % 形成等效动力学刚度矩阵
    K_eff = K + a0 * M + a1 * C;

    % ------------------- 关键修正开始 -------------------
    % 正确进行 LU 分解 (带行置换矩阵 P)
    [L, U, P] = lu(K_eff);   % 此时 P * K_eff = L * U
    % --------------------------------------------------

    % 时间步进循环
    for i = 1:(Nt - 1)
        % 计算 t(i+1) 时刻的等效载荷向量
        F_eff = F(:, i+1) ...
                + M * (a0 * X(:,i) + a2 * V(:,i) + a3 * A(:,i)) ...
                + C * (a1 * X(:,i) + a4 * V(:,i) + a5 * A(:,i));

        % 求解位移: 先置换, 再前代, 最后回代
        % 解 P * K_eff * X = P * F_eff  ->  (L*U) * X = P * F_eff
        X(:, i+1) = U \ (L \ (P * F_eff));

        % 更新加速度和速度 (保持原有 Newmark 公式)
        A(:, i+1) = a0 * (X(:, i+1) - X(:,i)) - a2 * V(:,i) - a3 * A(:,i);
        V(:, i+1) = V(:,i) + dt * ((1 - gamma) * A(:,i) + gamma * A(:, i+1));
    end

    % 输出 (与输入 t 相同)
    t = t(:)';   % 确保行向量
end