clc;clear
%% 和低压电对比频率
No_INTpoint_x=4;
No_INTpoint_y=4;
No_INTpoint_z=2;

% Alu plate（各向同性）
densityb=2700;
E_b = 70e9; nu_b = 0.33; H_b = 0.2e-3;

densityp=5440;
% 压电材料性材料每行为[E1, E2, G, nu12, nu21，h]
middle_pian = [
    30.336e9, 15.857e9, 5.515e9, 0.31, 0.162 ,0.3e-3];
e31=-6.5;e32=-6.5;e24=0;e15=0;E33=1.3e-8;
d31=-1.7e-10;d33=4e-10;     %(1000vpp左右)
d31=-1.4e-10; d33=3e-10;    %(500vpp)
d31=-1.3e-10; d33=2.51e-10; %(200vpp)

% Epoxy and pizplate（各向同性）
densitye=1100;
E_e = 3e9; nu_e = 0.32;
G_e=E_e/(2*(1+nu_e));

%damping
zeta=0.01378;
eta_b = 2*zeta;  % 铝板的损耗因子
eta_p = 2*zeta;  % 压电材料的损耗因子 

%流体
rho_fluid=1000;
mu = 1.002e-3; % 水的动力粘度 (Pa.s)


% a11=-3.9429;
% a22=100.0;
% b11= -17.6199;
% a11=0.1957;
% a22=20;
% b11= 20;

%厚度比
hp=0.3e-3;
hb=0.2e-3; 
he = 0.15e-3;
H=2*hp+hb+2*he;

bt=1;                        at=(0.5*hb+he)/(0.5*H);
b1=(0.5*hb+he)/(0.5*H);     a1=(0.5*hb)/(0.5*H);  
bm=(0.5*hb)/(0.5*H);         am=-(0.5*hb)/(0.5*H);
b2= -(0.5*hb)/(0.5*H);    a2=-(0.5*hb+he)/(0.5*H);
bb=-(0.5*hb+he)/(0.5*H);        ab=-1;

bp=[bt,bb,bt,bb];
ap=[at,ab,at,ab];

%% 单元编号
[nodes, elements,elements_index] = read_comsol8('5.16wangge.bdf');
nodes_ext=1e-3*nodes;
combined_index = [elements_index{2}];
elements_p = elements(combined_index, :);

%% 节点的物理坐标
% 上层节点坐标
node_total = size(nodes_ext,1);  %总节点数
coop = zeros(node_total, 3);     %初始化坐标矩阵，避免索引越界
coop(:,1:2)=nodes_ext(:,1:2);
coop(:,3)=0.5*H;

% 下层节点坐标
cop = zeros(node_total, 3);
cop(:,1:2)=nodes_ext(:,1:2);
cop(:,3)=-0.5*H;

disp_node = zeros(node_total, 5);  %初始化节点自由度矩阵
disp_node(1:node_total, 1:5) = 1;
% 固定部分节点自由度为0（约束）
zuobiaox = find(abs(nodes_ext(:,1)) < 1e-6);
disp_node(zuobiaox, :) = 0;

dof = 0;  %总自由度计数
for ni=1:node_total
    for nj=1:5
        if disp_node(ni,nj) ~= 0
            dof = dof + 1;
            disp_node(ni,nj) = dof;  %分配自由度编号
        end
    end
end

%节点平均坐标
cooo = (coop + cop) / 2;
x1 = coop(:,1); y1 = coop(:,2); z1 = coop(:,3);
x2 = cop(:,1); y2 = cop(:,2); z2 = cop(:,3);
%% 压电片的数量与位置
np=2;
%nn=elements_index{1};
nn=combined_index;
sys_dof = dof;
k = zeros(sys_dof, sys_dof);    %系统刚度矩阵 
K_complex = zeros(sys_dof, sys_dof); % 系统复刚度矩阵 (用于数值积分)
m = zeros(sys_dof, sys_dof);    %系统质量矩阵 
ka = zeros(sys_dof, np); 
ks = zeros(np, sys_dof); 

%% 结构刚度、质量矩阵组装 & 缓存水动力所需变量
jdzbp = coop;
jdzb1p = cop;
dybh = elements;

% 预分配元胞数组，用于给后面的水动力组装缓存必要数据
zmtemp_cell = cell(size(elements,1), 1); % 缓存坐标系矩阵
id_cell = cell(size(elements,1), 1);     % 缓存单元对应的全局非零索引
mask_cell = cell(size(elements,1), 1);   % 缓存单元非零索引的掩码

for loopi = 1:size(elements,1)
    isMFC = ismember(loopi, nn);
    % --- 1. 结构基板部分 ---
    [ek_b, ~, xv2i, xv1i, ~, zmtemp, v3i, ~, jtemp] = middleshellek(am, bm, E_b, nu_b, H, loopi, jdzbp, jdzb1p, dybh);
    zmtemp_cell{loopi} = zmtemp; 
    [em_b] = middleshellem(am, bm, zmtemp, v3i, densityb, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
    
    ek_total = ek_b;
    em_total = em_b;
    ek_complex_local = ek_b * (1 + 1i * eta_b);
    
    % --- 2. MFC压电片部分 ---
    if isMFC
        [ek_m1, D] = middleshellek_p(at, bt, middle_pian, H, loopi, jdzbp, jdzb1p, dybh);
        [em_m1] = middleshellem(at, bt, zmtemp, v3i, densityp, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        [ek_m2] = middleshellek_p(ab, bb, middle_pian, H, loopi, jdzbp, jdzb1p, dybh);
        [em_m2] = middleshellem(ab, bb, zmtemp, v3i, densityp, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        
        [ek_e1, ~] = middleshellek(a1, b1, E_e, nu_e, H, loopi, jdzbp, jdzb1p, dybh);
        [em_e1] = middleshellem(a1, b1, zmtemp, v3i, densitye, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        [ek_e2, ~] = middleshellek(a2, b2, E_e, nu_e, H, loopi, jdzbp, jdzb1p, dybh);
        [em_e2] = middleshellem(a2, b2, zmtemp, v3i, densitye, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        
        ek_total = ek_total + ek_m1 + ek_m2 + ek_e1 + ek_e2;
        em_total = em_total + em_m1 + em_m2 + em_e1 + em_e2;
        % 阻尼（MFC区域）
        ek_complex_local = ek_total * (1 + 1i * eta_p);
    end
    
    % --- 3. 提取自由度索引并缓存 ---
    index = zeros(40, 1);
    for zi = 1:8
        index((zi-1)*5 + 1) = disp_node(elements(loopi, zi), 1);
        index((zi-1)*5 + 2) = disp_node(elements(loopi, zi), 2);
        index((zi-1)*5 + 3) = disp_node(elements(loopi, zi), 3);
        index((zi-1)*5 + 4) = disp_node(elements(loopi, zi), 4);
        index((zi-1)*5 + 5) = disp_node(elements(loopi, zi), 5);
    end
    nonzero_idx = index(index ~= 0);
    valid_mask = index ~= 0;
    
    % 将索引和掩码存入缓存，迭代湿模态时无需再算这步
    id_cell{loopi} = nonzero_idx;
    mask_cell{loopi} = valid_mask;
    
    % --- 4. 组装全局结构刚度和质量矩阵 ---
    k(nonzero_idx, nonzero_idx) = k(nonzero_idx, nonzero_idx) + ek_total(valid_mask, valid_mask);
    K_complex(nonzero_idx, nonzero_idx) = K_complex(nonzero_idx, nonzero_idx) + ek_complex_local(valid_mask, valid_mask);
    m(nonzero_idx, nonzero_idx) = m(nonzero_idx, nonzero_idx) + em_total(valid_mask, valid_mask);
end

%% 组装压电耦合矩阵ka与ks
nnd=[nn,nn];
for n=1:np
    for loopi=1:size(elements,1)
        EE = ismember(loopi, nnd(:,n));
        if EE ~= 0
            dyhm = loopi;
            [ekat_M] = MFC_topka(ap(1,n), bp(1,n), H, D, hp, d33, d31, dyhm, jdzbp, jdzb1p, dybh);
            eka = ekat_M;
        else
            eka = zeros(40, 1);
        end
        for zi=1:8
            index((zi-1)*5 + 1) = disp_node(elements(loopi, zi), 1);
            index((zi-1)*5 + 2) = disp_node(elements(loopi, zi), 2);
            index((zi-1)*5 + 3) = disp_node(elements(loopi, zi), 3);
            index((zi-1)*5 + 4) = disp_node(elements(loopi, zi), 4);
            index((zi-1)*5 + 5) = disp_node(elements(loopi, zi), 5);
        end
        % 组装耦合矩阵
        for jx=1:40
            if index(jx) ~= 0
                ka(index(jx), n) = ka(index(jx), n) + eka(jx, 1);       
            end
        end
    end
end
%% 模态分析(空气)
nm=20;
[v,d]=eigs(k,m,nm,'SM');
tempd = diag(d);
[d, sortindex] = sort(tempd);
v = v(:, sortindex);
omega=(sqrt(d));
frequency = (sqrt(d))/(2*pi);  %转换频率(Hz)


%模态空间
nn=1;%截断模态
v=v(:,nn);
% Mw_modal = v' * Mw_global * v;
% Cw_modal = v' * Cw_global * v;
% mi1_global= v' * mi1_global * v;
K_modal  = v' * k * v;
K_complex_modal  = v' *  K_complex * v;
M_modal  = v' * m * v;

%% --- 实验标定参数输入区 ---
% 1. 仅使用刚度比例阻尼 (根据自由衰减实验的一阶衰减结果填入)
f_ex = frequency(1);  
omega_ex = 2 * pi * f_ex;              % 圆频率
beta_rayleigh = 2 * zeta / omega_ex;
C_matrix1 = beta_rayleigh * k;         % 真实的粘性阻尼矩阵

C_matrix2 = imag(K_complex) / omega_ex; %复模态阻尼
C_modal=v' * C_matrix2  * v;

% 2. 压电非线性饱和与迟滞等效参数
alpha_sat = 1;                   % 300V下的极化饱和衰减系数 (建议调试范围: 0.3 ~ 0.5)
psi_hys_deg = 15;                % 迟滞环导致的力学滞后角 (建议调试范围: 10° ~ 20°)
psi_hys = psi_hys_deg * (pi / 180); % 转为弧度


% 3. 修正制动厚度
h_act = 0.18e-3;                    % 实际真实的压电层厚度 (mm)
gamma_geo = h_act / hp ;
ka_true = ka * gamma_geo;

% 4.修正机电耦合矩阵 (降输入)
ka_eff = ka_true * alpha_sat;

%% --- 时间历程与非线性激振力生成区 ---
f_ex = frequency(1);             % 锁定一阶模态频率作为激励频率
V0 = 100;                       
Va = [V0; -V0];                    % 交流电压幅值向量 [压电片1; 压电片2]
Va_D = [0; 0];                   % 直流偏置电压向量 [压电片1; 压电片2]
phi = [0; 0];                    % 相位角差向量 [压电片1; 压电片2]

% 设定时间步和总时长
cycles = 100;                    % 仿真 100 个周期以确保达到稳态
T_end = cycles * (1 / f_ex);
points_per_cycle = 40;           % 每个周期采 40 个点保证精度
dt = (1 / f_ex) / points_per_cycle;
t = 0:dt:T_end;                  % 时间向量 (尺寸: 1 x Nt)


% 在交流电压中引入物理滞后角 psi_hys
V_ac_eff = Va .* sin(phi + 2 * pi * f_ex * t ); 
V_total_eff =  V_ac_eff+Va_D;

% 生成修正后的真实激振力矩阵
F_matrix1 = ka_eff * V_total_eff; %[系统自由度, 时间步数]     
F_matrix2 = ka_eff*(Va+Va_D);     %[系统自由度, 激励数量]   


F0=F_matrix2;
X_static =k \F0;
[DIS_D,dip_D]= ttq(disp_node, X_static, nodes_ext);
figure(11)
faces4 = elements(:, [1 2 3 4]);
patch('Faces', faces4, 'Vertices', DIS_D, ...
'FaceVertexCData', abs(dip_D), ...
'FaceColor', 'interp', 'EdgeColor', 'none');
axis equal; colorbar; colormap(turbo)
title('挠曲变形图')


%% 优化
f_exp = [5, 6, 7, 8, 9, 10]; % 实验频率点 (Hz)
%amp_exp = [2.153, 2.964, 3.004, 2.537, 1.978, 1.567] * 1e-3; % 300v双片对应的稳态振幅 (m)
amp_exp = [0.568, 0.984, 1.399, 1.209, 0.802, 0.498] * 1e-3; %  100v双片对应的稳态振幅 (m)
jth=8;

% 初始猜测值
initial_guess = [-5, -10, -5];
% 设定上下界 
lb = [-50, -50, -50]; 
ub = [ 10, 500,  100];

% 边界也只留两个！彻底放开 a2 的下界
params_init = [-5, -5]; 
lb = [-30, -30]; 
ub = [ 10,  10];



% 定义匿名目标函数 (固定其他已知参数)
obj_func = @(p) fit_fluid_params2( p, f_exp, amp_exp,  m, F_matrix2, jth, elements, zmtemp_cell, rho_fluid, mu, disp_node, v, nn ,  M_modal, K_complex_modal);
% 设置优化选项 (显示迭代过程)
% options = optimoptions('patternsearch', 'Display', 'iter');
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'StepTolerance', 1e-10);
% 运行拟合！
fprintf('开始参数辨识...\n');
fprintf('启动无导数算法进行参数辨识...\n');
% best_params = patternsearch(obj_func, initial_guess, [], [], [], [], lb, ub, [], options);
best_params = fmincon(obj_func, initial_guess, [], [], [], [], lb, ub, [], options);

% 提取结果
% opt_a1 = best_params(1);
% opt_a2 = best_params(2);
% opt_b1 = best_params(3);

opt_a1 = best_params(1);
opt_a2 = 0;
opt_b1 = best_params(2);
fprintf('拟合成功！\n最优 a11: %.4f\n最优 a22: %.4f\n最优 b11: %.4f\n', opt_a1, opt_a2, opt_b1);


%% --- 1. 结构湿模态求解  ---
tol = 1e-4; % 迭代误差容限
max_iter = 50; % 最大迭代次数
% 1. 给一个初始猜测频率（比如用你的实验共振频率，或者干模态频率）
f_guess = 7.0; % Hz
omega_guess = 2 * pi * f_guess;

%2. 水动力矩阵组装
fprintf('开始迭代求解真实湿模态频率...\n');
for iter = 1:max_iter
    
    % 每次迭代前，必须先将水动力矩阵清零
    Mw_global  = zeros(size(k, 1));
    Cw_global  = zeros(size(k, 1));
    mi1_global = zeros(size(k, 1));

    for loopi = 1:size(elements,1)
        % 直接从前面的一次性缓存中提取，省略重复计算
        zmtemp_local = zmtemp_cell{loopi};
        id           = id_cell{loopi};
        mask         = mask_cell{loopi};

        % 基于当前迭代步骤的 omega 计算单元水动力矩阵
        [mw_eff, cw_eff, mi1] = wateradding(zmtemp_local, rho_fluid, mu, omega_guess, opt_a1, opt_a2, opt_b1);

        % 直接组装至全局水动力矩阵
        Mw_global(id, id)  = Mw_global(id, id)  + mw_eff(mask, mask);
        Cw_global(id, id)  = Cw_global(id, id)  + cw_eff(mask, mask);
        mi1_global(id, id) = mi1_global(id, id) + mi1(mask, mask);
    end



    % 3. 组装总质量矩阵并求解新的特征值
    MM = m + Mw_global;
    [~, d] = eigs(k, MM, 1, 'SM'); % 只求一阶
    omega_new = sqrt(real(d(1,1)));
    f_new = omega_new / (2*pi);

    % 4. 检查是否收敛
    error_f = abs(f_new - f_guess);
    fprintf('Iter %d: 猜测频率 = %.4f Hz, 算出频率 = %.4f Hz, 误差 = %.6f\n', ...
        iter, f_guess, f_new, error_f);

    if error_f < tol
        fprintf('求解成功！真实的理论水下谐振频率为: %.4f Hz\n', f_new);
        break;
    end

    % 5. 更新猜测值，继续迭代
    f_guess = f_new;
    omega_guess = omega_new;
end



% %% 1. 载入你刚得到的最优参数(频域)
% a1_opt = -3.3745;
% a2_opt = 0;       % 已降阶，严格锁定为 0
% b1_opt = -6.2271;
% 
% %% 2. 迭代计算包含“附加质量”的真实水下共振频率
% % (这段逻辑你在前面已经跑通了，这里直接用你最后收敛的频率即可)
% f_wet = f_new; % 你的水下谐振频率 (Hz)
% omega_wet = 2 * pi * f_wet;
% 
% %% 3. 在该共振频率下，组装水动力矩阵 
% Mw_global = zeros(sys_dof);
% Cw_global = zeros(sys_dof);
% 
% for loopi = 1:size(elements,1)
%     zmtemp_local = zmtemp_cell{loopi};
%     id   = id_cell{loopi};
%     mask = mask_cell{loopi};
% 
%     % 调用你的水动力函数 (传入最优参数)
%     [mw_eff, cw_eff] = wateradding(zmtemp_local, rho_fluid, mu, omega_wet, a1_opt, a2_opt, b1_opt);
% 
%     Mw_global(id, id) = Mw_global(id, id) + mw_eff(mask, mask);
%     Cw_global(id, id) = Cw_global(id, id) + cw_eff(mask, mask);
% end
% 
% % 投影到模态空间
% Mw_modal = v' * Mw_global * v;
% Cw_modal = v' * Cw_global * v;
% 
% % 计算该频率下的全局动态刚度矩阵 Z
% Z_wet = K_modal - (omega_wet^2) * (M_modal + Mw_modal) + 1i * omega_wet * (C_modal + Cw_modal);
% 
% %% 4. 模拟 100V - 500V 的稳态时域响应
% V_list = [100, 200, 300, 400, 500];
% t = linspace(0, 1, 1000); % 模拟 1 秒钟的时域波形
% 
% figure('Color', 'w', 'Position', [100, 100, 800, 500]);
% hold on; grid on;
% colors = lines(length(V_list)); % 取 5 种不同的颜色
% 
% % 假设你有一个单位电压(1V)对应的力向量 F_amp_unit
% F_amp_unit = F_matrix2 / 100; 
% 
% sim_amps = zeros(1, length(V_list)); % 用于记录各电压下的峰值振幅
% 
% for i = 1:length(V_list)
%     V_current = V_list(i);
% 
%     % 1. 计算当前电压下的模态力
%     F_modal_current =v' * (F_amp_unit * V_current);
% 
%     % 2. 求解复数频域响应 (稳态)
%     q_amp = Z_wet \ F_modal_current;
% 
%     % 3. 映射回物理空间，提取尖端自由度的复数振幅
%     U_amp = v * q_amp;
%     tip_complex_amp = U_amp(jth) * 1000; % 转换为 mm
% 
%     % 4. 提取幅值和相位
%     amp_mm = abs(tip_complex_amp);
%     phase_rad = angle(tip_complex_amp);
%     sim_amps(i) = amp_mm;
% 
%     % 5. 重构时域稳态波形: w(t) = A * sin(wt + phi)
%     % (假设外力为 sin(wt)，如果是 cos，只需将下面的 sin 改为 cos 即可)
%     tip_time_domain = amp_mm * sin(omega_wet * t + phase_rad);
% 
%     % 6. 绘图
%     plot(t, tip_time_domain, 'LineWidth', 1.5, 'Color', colors(i,:), ...
%          'DisplayName', sprintf('%d V (Amp: %.2f mm)', V_current, amp_mm));
% end
% 
% xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Tip Displacement (mm)', 'FontSize', 12, 'FontWeight', 'bold');
% title(sprintf('Steady-State Time Response at %.2f Hz', f_wet), 'FontSize', 14);
% legend('Location', 'northeastoutside', 'FontSize', 11);
% xlim([0, 0.5]); % 只看前 0.5 秒的波形，方便看清正弦波
% set(gca, 'FontSize', 11);
% 
% fprintf('\n==== 各电压下模拟共振振幅 ====\n');
% for i = 1:length(V_list)
%     fprintf('%d V: %.3f mm\n', V_list(i), sim_amps(i));
% end






%% 1. 载入最优流体参数 (时域)
a1_opt = -3.3730;
a2_opt = 0;       
b1_opt = -13.3615;

f_wet = 7.3553; % 你的理论水下共振频率 (Hz)
omega_wet = 2 * pi * f_wet;

%% 2. 组装物理水动力矩阵
Mw_global = zeros(sys_dof);
Cw_global = zeros(sys_dof);

for loopi = 1:size(elements,1)
    zmtemp_local = zmtemp_cell{loopi};
    id   = id_cell{loopi};
    mask = mask_cell{loopi};
    [mw_eff, cw_eff] = wateradding(zmtemp_local, rho_fluid, mu, omega_wet, a1_opt, a2_opt, b1_opt);
    Mw_global(id, id) = Mw_global(id, id) + mw_eff(mask, mask);
    Cw_global(id, id) = Cw_global(id, id) + cw_eff(mask, mask);
end

%% 3. 将所有矩阵投影到模态空间
% 将物理矩阵转换为降阶后的模态矩阵，极大缩短 Newmark 积分时间！
MM_modal = M_modal + v' * Mw_global * v;
C_matrix2 = imag(K_complex) / omega_wet; %复模态阻尼
C_modal=v' * C_matrix2  * v;
CC_modal = C_modal + v' * Cw_global * v;
KK_modal = K_modal;

%% 4. 时域积分设置
dt = 0.001;  % 时间步长 1ms
t_end = 4.0; % 模拟 4 秒，确保进入稳态
t = 0:dt:t_end;

% 提取单位电压对应的“模态力”向量
% (假设 F_amp_vector 是物理空间 100V 对应的力向量)
F_amp_unit = F_matrix2 / 100; 
F_modal_unit = v' * F_amp_unit; 

V_list = [100, 200, 300, 400, 500];
steady_amps = zeros(1, length(V_list));

figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on; grid on;
colors = lines(length(V_list));

%% 5. 循环计算 100V 到 500V 的时域响应
for v_idx = 1:length(V_list)
    V_current = V_list(v_idx);
    
    % 构造时间历程激振力矩阵 (尺寸: [模态数, 时间步数])
    % 这里是纯正弦波激振力 (后续如果加上迟滞非线性，改这里即可)
    F_matrix = (F_modal_unit * V_current) * sin(omega_wet * t); 
    
    % 调用你自己的积分器！（在模态空间内极速求解）
    fprintf('\n开始 Newmark-Beta 时域积分 (当前激励: %d V，共 %d 个时间步)...\n', V_current, length(t));
    [t_out1, q_out] = Standard_Newmark_Beta2(MM_modal, CC_modal, KK_modal, F_matrix, t);
    fprintf('积分完成！\n');
    
    % 将计算出的模态坐标 q_out 映射回物理空间的位移 X_out
    X_out = v * q_out; 
    
    % 提取尖端自由度的位移，并转换为 mm
    tip_disp_mm = X_out(jth, :) * 1000;
    
    % 提取后 30% 时长的数据来计算“稳态振幅”
    steady_part = tip_disp_mm(floor(0.7 * length(t)) : end);
    steady_amps(v_idx) = (max(steady_part) - min(steady_part)) / 2;
    
    % 绘制该电压下的时域波形
    plot(t_out1, tip_disp_mm, 'LineWidth', 1.2, 'Color', colors(v_idx,:), ...
         'DisplayName', sprintf('%d V (稳态: %.2f mm)', V_current, steady_amps(v_idx)));
end

%% 6. 绘图修饰
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Tip Displacement (mm)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('True Time-Domain Response at %.2f Hz', f_wet), 'FontSize', 14);
legend('Location', 'northeastoutside', 'FontSize', 11);
set(gca, 'FontSize', 11);
xlim([0, 4.0]); 

fprintf('\n==== 基于纯时域积分 (Newmark-Beta) 的稳态振幅预估 ====\n');
for i = 1:length(V_list)
    fprintf('%d V: %.3f mm\n', V_list(i), steady_amps(i));
end












q_static = KK_modal \ (v' * F_matrix2 / 100 * 100); % 假设 100V
disp(['静态模态位移 q_static: ', num2str(q_static)]);

X_modal_to_physical = v * q_static; 
tip_displacement_from_modal = X_modal_to_physical(jth); % 模态算出来的尖端位移
tip_displacement_from_physical = X_static(jth);       % 物理空间直接算出来的尖端位移

fprintf('物理空间直接算出的尖端位移: %e\n', tip_displacement_from_physical);
fprintf('从模态坐标反推的尖端位移: %e\n', tip_displacement_from_modal);
fprintf('二者比值 (Scaling Factor): %e\n', tip_displacement_from_physical / tip_displacement_from_modal);

MM = m + Mw_global;           % 总质量
CC = C_matrix2 + Cw_global;    % 总阻尼
[t_out1, X_out] = Standard_Newmark_Beta(MM, CC, k, F_matrix1, t);

% 绘制特定自由度的时域响应 
jth=8; % 观察的自由度8
time_vector = 0:dt:T_end;
figure(3)                   
plot(t ,X_out(jth,:)*1e3)
xlabel('Time (seconds)')
ylabel('displacement(mm)')
hold on
% ylim([-2.5, 2.5]);
xlim([min(time_vector), max(time_vector)]);

time_vector= 2.5:dt:3.0;
figure(4)                   
plot(time_vector,X_out(jth,2.5/dt:3.0/dt)*1e3);%401:801
xlabel('Time (seconds)')
ylabel('displacement(mm)')
hold on
xlim([min(time_vector), max(time_vector)]);