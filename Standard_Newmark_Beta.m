function [t, X] = Standard_Newmark_Beta(M, C, K, F, t)
    % Standard_Newmark_Beta: 实数域 Newmark-beta 时域积分
    % 输入:
    %   M, C, K - 系统的质量、阻尼、实数刚度矩阵
    %   F - 激振力矩阵 (每一列对应一个时间步的力向量)
    %   t - 时间向量
    % 输出:
    %   t - 时间向量
    %   X - 位移响应矩阵 (每一列对应一个时间步的位移向量)
    
    % Newmark-beta 算法参数 (平均加速度法，无条件稳定)
    gamma = 0.5;
    beta = 0.25;
    
    dt = t(2) - t(1);
    Nt = length(t);
    N_dof = size(M, 1);
    
    % 初始化响应矩阵
    X = zeros(N_dof, Nt);    % 位移
    V = zeros(N_dof, Nt);    % 速度
    A = zeros(N_dof, Nt);    % 加速度
    
    % 初始条件 (假设静止起振)
    X(:,1) = zeros(N_dof, 1);
    V(:,1) = zeros(N_dof, 1);
    
    % 计算 t=0 时刻的初始加速度
    % M*A + C*V + K*X = F  =>  A = M \ (F - C*V - K*X)
    A(:,1) = M \ (F(:,1) - C*V(:,1) - K*X(:,1));
    
    % 计算等效动力学刚度矩阵 (仅需计算一次)
    a0 = 1 / (beta * dt^2);
    a1 = gamma / (beta * dt);
    K_eff = K + a0 * M + a1 * C;
    
    % 对等效刚度矩阵进行一次 LU 分解，极大地提升循环内的求解速度
    [L_K, U_K] = lu(K_eff);
    
    % 积分常数预计算
    a2 = 1 / (beta * dt);
    a3 = 1 / (2 * beta) - 1;
    a4 = gamma / beta - 1;
    a5 = (dt / 2) * (gamma / beta - 2);
    
    % 时间步进循环
    for i = 1:(Nt - 1)
        % 计算 t(i+1) 时刻的等效载荷向量
        F_eff = F(:, i+1) + M * (a0 * X(:,i) + a2 * V(:,i) + a3 * A(:,i)) ...
                          + C * (a1 * X(:,i) + a4 * V(:,i) + a5 * A(:,i));
                      
        % 求解 t(i+1) 时刻的位移 (利用已有的 LU 分解)
        X(:, i+1) = U_K \ (L_K \ F_eff);
        
        % 更新 t(i+1) 时刻的加速度和速度
        A(:, i+1) = a0 * (X(:, i+1) - X(:,i)) - a2 * V(:,i) - a3 * A(:,i);
        V(:, i+1) = V(:,i) + dt * ((1 - gamma) * A(:,i) + gamma * A(:, i+1));
    end
    
    t_out = t;
end