function error_sum = fit_fluid_params2(params, f_exp, amp_exp, M_struct, F_amp_vector, tip_dof, ...
    elements, zmtemp_cell, rho_fluid, mu, disp_node, ...
    Phi,nn, M_modal, K_complex_modal)

% a1 = params(1);
% a2 = params(2);
% b1 = params(3);

a2 =0; 
a1 =params(1);
b1 = params(2);

error_sum = 0;
Phi=Phi(:,nn);


sys_dof = size(M_struct,1);

% ===========================
% modal force projection
% ===========================
F_modal = Phi' * F_amp_vector;

for k_idx = 1:length(f_exp)

    omega = 2*pi*f_exp(k_idx);

    Mw_global = zeros(sys_dof);
    Cw_global = zeros(sys_dof);

    % ======================================================
    % 1. assemble fluid operator in physical space
    % ======================================================
    for loopi = 1:size(elements,1)

        zmtemp_local = zmtemp_cell{loopi};

        [mw_eff, cw_eff] = wateradding( ...
            zmtemp_local, rho_fluid, mu, omega, a1, a2, b1);

        index = zeros(40,1);

        for zi=1:8
            index((zi-1)*5+1) = disp_node(elements(loopi,zi),1);
            index((zi-1)*5+2) = disp_node(elements(loopi,zi),2);
            index((zi-1)*5+3) = disp_node(elements(loopi,zi),3);
            index((zi-1)*5+4) = disp_node(elements(loopi,zi),4);
            index((zi-1)*5+5) = disp_node(elements(loopi,zi),5);
        end

        id = index(index>0);
        mask = index~=0;

        Mw_global(id,id) = Mw_global(id,id) + mw_eff(mask,mask);
        Cw_global(id,id) = Cw_global(id,id) + cw_eff(mask,mask);

    end

    % ======================================================
    % 2. PROJECT TO MODAL SPACE
    % ======================================================
    Mw_modal = Phi' * Mw_global * Phi;
    Cw_modal = Phi' * Cw_global * Phi;
    if real(Cw_modal) <= 0
        % 如果计算出的附加水阻尼小于等于 0，触发惩罚机制！
        % 给一个极大的误差值 (例如 1e10)，警告优化器此路不通
        error_sum = 1e10;
        % 直接退出当前函数，不需要再做后续无意义的矩阵求解，极大节省算力！
        return;
    end
    % ======================================================
    % 3. MODAL COUPLED DYNAMIC STIFFNESS
    % ======================================================
    Z = K_complex_modal ...
        - (omega^2) * (M_modal + Mw_modal) ...
        + 1i * omega * Cw_modal;

    % ======================================================
    % 4. SOLVE IN MODAL SPACE
    % ======================================================
    q = Z \ F_modal;

    % back to physical space
    U = Phi * q;

    sim_amp = abs(U(tip_dof));


    % ======================================================
    % 5. error
    % ======================================================
    error_sum = error_sum + ((sim_amp - amp_exp(k_idx))*1000)^2;

end
end