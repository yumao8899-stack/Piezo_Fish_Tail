function error_sum = fit_fluid_params(params, f_exp, amp_exp, M_struct, K_struct, C_matrix, F_amp_vector, tip_dof, elements, zmtemp_cell, rho_fluid, mu, disp_node)
    a1 = params(1);
    a2 = params(2);
    b1 = params(3);
    
    error_sum = 0;
    sys_dof = size(M_struct, 1);
    
    for k_idx = 1:length(f_exp)
        omega = 2 * pi * f_exp(k_idx);
        
        Mw_global = zeros(sys_dof, sys_dof);
        Cw_global = zeros(sys_dof, sys_dof);
        
        for loopi = 1:size(elements, 1)
            zmtemp_local = zmtemp_cell{loopi}; 
            
            % 调用水动力基底矩阵
            [mw_eff, cw_eff] = wateradding(zmtemp_local, rho_fluid, mu, omega, a1, a2, b1);
            
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
        
        % -----------------------------------------------------------------
        % 最终系统矩阵组装 (基于纯实数刚度和 C_matrix)
        % -----------------------------------------------------------------
        M_total = M_struct + Mw_global;
        C_total = C_matrix + Cw_global; % 结构瑞利阻尼与流体阻尼直接叠加
        
        % 标准频域代数方程: (-w^2 M + iw C + K) * U = F
        H_dyn = K_struct - (omega^2) * M_total + (1i * omega) * C_total;
        
        % 求解稳态位移
        U_complex = H_dyn \ F_amp_vector;
        
        % 提取靶点振幅
        sim_amp = abs(U_complex(tip_dof));
        
        % L2 范数累加
        error_sum = error_sum + ((sim_amp - amp_exp(k_idx)) * 1000)^2;
    end
end