function [f_vec, X_complex, X_amp, X_phase] = Harmonic_Response_Solver(M, K_complex, F_amp, f_vec)
    % Harmonic_Response_Solver: 结构频域（谐响应）求解器
    % 
    % 输入参数:
    %   M         - 系统的全局质量矩阵 (N_dof x N_dof)
    %   K_complex - 系统的全局复刚度矩阵 (N_dof x N_dof)，已包含阻尼 K * (1 + 1i * eta)
    %   F_amp     - 激振力的复幅值向量 (N_dof x 1)，表示稳态激励力
    %   f_vec     - 需要扫频的频率向量 (1 x Nf) (单位: Hz)
    %
    % 输出参数:
    %   f_vec     - 频率向量 (与输入相同)
    %   X_complex - 复数位移响应矩阵 (N_dof x Nf)
    %   X_amp     - 位移幅值矩阵 (N_dof x Nf)
    %   X_phase   - 位移相位矩阵 (N_dof x Nf，单位: 弧度)
    
    % 获取系统自由度和频率点数
    N_dof = size(M, 1);
    Nf = length(f_vec);
    
    % 预分配内存以提升计算速度
    X_complex = zeros(N_dof, Nf);
    
    % 遍历每一个计算频率点
    for i = 1:Nf
        omega = 2 * pi * f_vec(i); % 将频率转换为圆频率 (rad/s)
        
        % -------------------------------------------------------------
        % 核心动力学方程: (-omega^2 * M + K_complex) * X = F_amp
        % -------------------------------------------------------------
        % 构建动刚度矩阵 (Dynamic Stiffness Matrix)
        D_matrix = -omega^2 * M + K_complex; 
        
        % 求解代数方程组，得到该频率下的复数稳态响应幅值
        % 警告：不要用 inv(D_matrix)*F_amp，使用左除 (\) 速度更快且精度更高
        X_complex(:, i) = D_matrix \ F_amp; 
    end
    
    % 从复数结果中提取物理意义直观的幅值和相位
    X_amp = abs(X_complex);
    X_phase = angle(X_complex);
end