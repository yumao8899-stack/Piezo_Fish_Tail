function [f, Hj] = calculate_frf(~, m, ka, Va, jth, omega, v, ac, bc, nm)
%扫频位移响应（基于模态叠加法）
% 输入:
%   k, m    - 刚度、质量矩阵
%   ka      - 压电耦合矩阵 (ndof x np)
%   Va      - 电压向量 (np x 1)，如 [600; -600]
%   jth     - 观测自由度索引
%   omega   - 模态圆频率 (从 eigs 得到)
%   v       - 模态矩阵
%   ac, bc  - Rayleigh 阻尼系数
%   nm      - 使用的模态数
% 输出:
%   f       - 扫频频率 [Hz]
%   U_amp   - 位移幅值 [m]（jth 自由度的稳态响应）

    %% 参数设置
    f = linspace(0, 1000, 5000);  % 扫频范围 0~1000 Hz
    w = 2*pi * f;
    nf = length(f);
    
    %% 质量归一化
    V_use = v(:, 1:nm);
    omega_use = omega(1:nm);
    Mmodal = V_use' * m * V_use;
    Vnorm = V_use * diag(1./sqrt(diag(Mmodal)));
    
    %% 模态阻尼
    zeta = 0.5 * (ac./omega_use + bc.*omega_use);
    %% 模态力
    F0 = ka * Va;
    Fr = Vnorm' * F0;  % (nm x 1)
    %% 扫频计算
    Hj = zeros(nf, 1);
    for i = 1:nf
        U_modal = 0;
        for r = 1:nm
            % 第 r 阶模态在频率 w(i) 处的贡献
            phi_jr = Vnorm(jth, r);
            % 模态响应传递函数
            H_r = 1 / (omega_use(r)^2 - w(i)^2 + 1i*2*zeta(r)*omega_use(r)*w(i));        
            % 累加到物理坐标
            U_modal = U_modal + phi_jr * Fr(r) * H_r;
        end
        Hj(i) = abs(U_modal);  % 取幅值
    end
end