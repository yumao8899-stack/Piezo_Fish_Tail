function [mw_eff, cw_eff,mi1] = wateradding(zmtemp, rho_fluid, mu, omega, a1, a2, b1)
   
%%zmtemp是中性面的各点坐标

% 仅在第3自由度 (Z向横向弯曲挠度 w) 上起作用的附加水质量矩阵
    mw_eff = zeros(40, 40);
    cw_eff = zeros(40, 40);
    mi1= zeros(40, 40);

    % 2x2 二维高斯积分点
    gauss_pt = [-0.774596669241483, 0.0, 0.774596669241483];
    gauss_wt = [0.555555555555556, 0.888888888888889, 0.555555555555556];
    
    for i = 1:3
        for j = 1:3
            zb1 = gauss_pt(i); 
            zb2 = gauss_pt(j);
            W   = gauss_wt(i) * gauss_wt(j); % 提取高斯权重
            
            % 形函数
            ni(1)=1/4*(1+zb1)*(1+zb2)*(zb1+zb2-1);
            ni(2)=1/4*(1-zb1)*(1+zb2)*(-zb1+zb2-1);
            ni(3)=1/4*(1-zb1)*(1-zb2)*(-zb1-zb2-1);
            ni(4)=1/4*(1+zb1)*(1-zb2)*(zb1-zb2-1);
            ni(5)=1/2*(1-zb1^2)*(1+zb2);
            ni(6)=1/2*(1-zb1)*(1-zb2^2);
            ni(7)=1/2*(1-zb1^2)*(1-zb2);
            ni(8)=1/2*(1+zb1)*(1-zb2^2);
             
            % 形函数各自偏导
            pni(1,1)=1/4*(1+zb2)*(zb1+zb2-1)+(1/4+1/4*zb1)*(1+zb2);
            pni(1,2)=(1/4+1/4*zb1)*(zb1+zb2-1)+(1/4+1/4*zb1)*(1+zb2);
            pni(2,1)=-1/4*(1+zb2)*(-zb1+zb2-1)-(1/4-1/4*zb1)*(1+zb2);
            pni(2,2)=(1/4-1/4*zb1)*(-zb1+zb2-1)+(1/4-1/4*zb1)*(1+zb2);
            pni(3,1)=-1/4*(1-zb2)*(-zb1-zb2-1)-(1/4-1/4*zb1)*(1-zb2);
            pni(3,2)=-(1/4-1/4*zb1)*(-zb1-zb2-1)-(1/4-1/4*zb1)*(1-zb2);
            pni(4,1)=1/4*(1-zb2)*(zb1-zb2-1)+(1/4+1/4*zb1)*(1-zb2);
            pni(4,2)=-(1/4+1/4*zb1)*(zb1-zb2-1)-(1/4+1/4*zb1)*(1-zb2);
            pni(5,1)=-zb1*(1+zb2);          pni(5,2)=1/2-1/2*zb1^2;
            pni(6,1)=-1/2+1/2*zb2^2;        pni(6,2)=-2*(1/2-1/2*zb1)*zb2;
            pni(7,1)=-zb1*(1-zb2);          pni(7,2)=-1/2+1/2*zb1^2;
            pni(8,1)=1/2-1/2*zb2^2;         pni(8,2)=-2*(1/2+1/2*zb1)*zb2;

            
            % 计算 2D 雅可比面积微元 dA 和该点的全局 X 坐标
            J11=0; J12=0; J21=0; J22=0; X_global=0;
            for k=1:8
                J11 = J11 + pni(k,1)*zmtemp(k,1);
                J12 = J12 + pni(k,1)*zmtemp(k,2);
                J21 = J21 + pni(k,2)*zmtemp(k,1);
                J22 = J22 + pni(k,2)*zmtemp(k,2);
                X_global = X_global + ni(k)*zmtemp(k,1);
            end
            dA = J11*J22 - J12*J21;
            
            % 1. 获取截面特征宽度
            b_local = get_fin_width(X_global);
           
            % 2. 计算特征频率 beta 
            beta_val = (rho_fluid * omega * b_local^2) / (2 * pi * mu);
            % 严密计算间隙宽弦比 delta (若没挖空，delta自然等于0)
            b_gap = get_gap_width(X_global);
            delta= b_gap / b_local;


            % 3. 计算实部与虚部 
            Re_Theta = 1.02 + a1*delta + (2.45 + a2*delta) * (beta_val^(-0.5));
            Im_mag   = (2.49 + b1*delta)* (beta_val^(-0.5)); % 虚部幅值
            

       
            % 质量面密度 = 原面密度 * Re[Theta]
            m_added_area = (pi/4) * rho_fluid * b_local * Re_Theta;
           
            % 阻尼面密度 = 原面密度 * Im_mag * omega
            % (因为阻尼对应力方程中的速度项 iw，提取出的 omega 必须乘在这里)
            c_added_area = (pi/4) * rho_fluid * b_local * Im_mag * omega;
            
            % 将水质量精准赋予第3自由度 (横向位移w)
            for node_m = 1:8
                for node_n = 1:8
                    idx_m = (node_m - 1) * 5 + 3;  
                    idx_n = (node_n - 1) * 5 + 3;
                    N_ij = ni(node_m) * ni(node_n) * dA * W;
                    mw_eff(idx_m, idx_n) = mw_eff(idx_m, idx_n) + m_added_area * N_ij;
                    cw_eff(idx_m, idx_n) = cw_eff(idx_m, idx_n) + c_added_area * N_ij;
                    mi1(idx_m, idx_n) = mi1(idx_m, idx_n) + (pi/4) * rho_fluid * b_local  * N_ij;
                end
            end
        end
    end
end



% =========================================================================
% 【子函数】：计算特征宽度 
% =========================================================================
function b_local = get_fin_width(x)
    bm = 0.049; be = 0.040; lm = 0.020; le = 0.055;
    
    if x <= bm
        b_local = lm;
    else
        % 归一化坐标 t，并添加边界保护
        t = min(max((x - bm) / be, 0), 1);
        % 立方形状函数 S(t) = 3t^2 - 2t^3
        St = 3 * t^2 - 2 * t^3;
        % b_local = lm + (le - lm) * S(t)
        b_local = lm + (le - lm) * St;
    end
end


% -------------------------------------------------------------
% 子函数 2：获取月牙间隙宽度 (Caudal Fork Gap)
% -------------------------------------------------------------
function b_gap = get_gap_width(x)
    L = 0.089; df = 0.010; le = 0.055;
    x_fork_start = L - df;
    
    if x <= x_fork_start
        b_gap = 0;
    else
        % 截断保护防止超出尾部
        x_eff = min(x, L);
        % 依据公式 x_tail(y) = L - df * [1 - (2y/le)^2]
        % 解得 gap = 2y = le * sqrt(1 - (L - x)/df)
        b_gap = le * sqrt(1 - (L - x_eff) / df);
    end
end
