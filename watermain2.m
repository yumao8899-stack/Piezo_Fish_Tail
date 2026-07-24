clc;clear
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
d31=-1.7e-10;d33=4e-10;%(1000vpp左右)
d31=-1.4e-10; d33=3e-10;%(500vpp)

%流体
rho_fluid=1000;
mu = 1.002e-3; % 水的动力粘度 (Pa.s)
omega=2*pi*7;
 % a1=1;  a2=1; b1=1;
% a1=1.54; a2=100; b1=-10.8;


a1= 1;
a2=100.0;
b1= -50.0;


eta_b = 0.0216;  % 铝板的损耗因子
eta_p = 0.0216;  % 压电材料的损耗因子 (若有独立测试数据可单独修改)

hp=0.3e-3;
hb=0.2e-3;
H=2*hp+hb;

bt=1;                      at=(0.5*hb)/(0.5*H);
bm=(0.5*hb)/(0.5*H);       am=-(0.5*hb)/(0.5*H);
bb=-(0.5*hb)/(0.5*H);      ab=-1;

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

Mw_global = zeros(sys_dof, sys_dof);
Cw_global = zeros(sys_dof, sys_dof);


ka = zeros(sys_dof, np); 
ks = zeros(np, sys_dof); 

%%
jdzbp = coop;
jdzb1p = cop;
dybh = elements;
index = zeros(40, 1);  %存储单元节点的系统自由度

zmtemp_cell = cell(size(elements,1), 1);%水动力元胞数组

% 组装刚度矩阵和质量矩阵
for loopi=1:size(elements,1)
    dyhm = loopi;
    BB = ismember(loopi, nn);
    [ekpp, ~, xv2i, xv1i, ~, zmtemp, v3i, ~, jtemp] = middleshellek(am, bm, E_b , nu_b, H, dyhm, jdzbp, jdzb1p, dybh);
    zmtemp_cell{loopi} = zmtemp;

    [empp] = middleshellem(am, bm, zmtemp, v3i, densityb, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
    ek_complex_local = ekpp * (1 + 1i * eta_b);
    if BB ~= 0
        [ekp1,D] = middleshellek_p(at,bt, middle_pian, H, dyhm, jdzbp, jdzb1p, dybh);
        [emm1] = middleshellem(at,bt, zmtemp, v3i, densityp, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        [ekp2] = middleshellek_p(ab,bb, middle_pian, H, dyhm, jdzbp, jdzb1p, dybh);
        [emm2] = middleshellem(ab,bb, zmtemp, v3i, densityp, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        ek =  ekp1 + ekp2 + ekpp;
        em =  emm1 + emm2 + empp;
        ek_complex_local = ek_complex_local + (ekp1 + ekp2) * (1 + 1i * eta_p);
    else
        ek =   ekpp;
        em =   empp;
    end
    [mw_eff, cw_eff] = wateradding(zmtemp, rho_fluid, mu, omega, a1, a2, b1);

    for zi=1:8
        index((zi-1)*5 + 1) = disp_node(elements(loopi, zi), 1);
        index((zi-1)*5 + 2) = disp_node(elements(loopi, zi), 2);
        index((zi-1)*5 + 3) = disp_node(elements(loopi, zi), 3);
        index((zi-1)*5 + 4) = disp_node(elements(loopi, zi), 4);
        index((zi-1)*5 + 5) = disp_node(elements(loopi, zi), 5);
    end
    % 组装整体矩阵
   % 3. 组装到全局矩阵
    nonzero_idx = index(index ~= 0);
    valid_mask = index ~= 0;
    
    % 向量化组装 
    k(nonzero_idx, nonzero_idx) = k(nonzero_idx, nonzero_idx) + ek(valid_mask, valid_mask);
    K_complex(nonzero_idx, nonzero_idx) = K_complex(nonzero_idx, nonzero_idx) + ek_complex_local(valid_mask, valid_mask);

    m(nonzero_idx, nonzero_idx) = m(nonzero_idx, nonzero_idx) + em(valid_mask, valid_mask);
    Mw_global(nonzero_idx, nonzero_idx) = Mw_global(nonzero_idx, nonzero_idx) + mw_eff(valid_mask, valid_mask);
    Cw_global(nonzero_idx, nonzero_idx) = Cw_global(nonzero_idx, nonzero_idx) + cw_eff(valid_mask, valid_mask);
end



%% 组装压电耦合矩阵ka与ks
nnd=[nn,nn];
for n=1:np
    for loopi=1:size(elements,1)
        EE = ismember(loopi, nnd(:,n));
        if EE ~= 0
            dyhm = loopi;
            [ekat_M] = MFC_topka(ap(1,n), bp(1,n), H, D, hp, d33, d31, dyhm, jdzbp, jdzb1p, dybh);
            eka = 2*ekat_M;
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

% (空气)
nm=20;
[v,d]=eigs(k,m,nm,'SM');
tempd = diag(d);
[d, sortindex] = sort(tempd);
v = v(:, sortindex);
omega=(sqrt(d));
frequency1 = (sqrt(d))/(2*pi);  %转换频率(Hz)


%模态空间
nn=1;%截断模态
v=v(:,nn);
Mw_modal = v' * Mw_global * v;
Cw_modal = v' * Cw_global * v;
K_modal  = v' * k * v;
M_modal  = v' * m * v;
%% --- 实验标定参数输入区 ---
f_ex = frequency1(1);  
% 1. 结构与流体耦合等效阻尼 (根据自由衰减实验的一阶衰减结果填入)
zeta_1 = 0.010588;                   % 假设实验测得的一阶阻尼比为 0.04
omega_ex = 2 * pi * f_ex;            % 圆频率

% 2. 压电非线性饱和与迟滞等效参数
alpha_sat = 0.3;                % 300V下的极化饱和衰减系数 (建议调试范围: 0.3 ~ 0.5)
psi_hys_deg = 15;               % 迟滞环导致的力学滞后角 (建议调试范围: 10° ~ 20°)
psi_hys = psi_hys_deg * (pi / 180); % 转为弧度

% 3. 修正制动厚度
h_act = 0.18e-3;                % 实际真实的压电层厚度 (mm)
gamma_geo = h_act / hp ;
ka_true = ka * gamma_geo;

V0 = -300;                       
Va = [V0; -V0];
%% --- 矩阵与激励重构区 ---
% 仅使用刚度比例阻尼 (假设单模态低频主导，忽略质量比例阻尼以防止刚体模态发散)
beta_rayleigh = 2 * zeta_1 / omega_ex;
C_matrix = beta_rayleigh * k; % 真实的粘性阻尼矩阵
C_modal=v' * C_matrix  * v;

% 修正机电耦合矩阵 (降输入)
ka_eff = ka_true * alpha_sat;
F_amp_vector = ka_eff * Va;
f_exp = [5, 6, 7, 8, 9, 10]; % 实验频率点 (Hz)
amp_exp = [2.153, 2.964, 3.004, 2.537, 1.978, 1.567] * 1e-3; % 对应的稳态振幅 (m)
jth=8;

%% 优化
% 初始猜测值 a1, a2, b1 通常在 -10 到 10 的范围内
initial_guess = [-5.0, -5.0, 1.0]; 

% 设定上下界 
lb = [-50, -100, -50]; 
ub = [ 50,  100,  50];

% 定义匿名目标函数 (固定其他已知参数)
obj_func = @(p) fit_fluid_params2( p, f_exp, amp_exp,  m, F_amp_vector, jth, elements, zmtemp_cell, rho_fluid, mu, disp_node, v, nn ,  M_modal, K_modal, C_modal);
% 设置优化选项 (显示迭代过程)
options = optimoptions('patternsearch', 'Display', 'iter');
% 运行拟合！
fprintf('开始参数辨识...\n');
fprintf('启动无导数算法进行参数辨识...\n');
best_params = patternsearch(obj_func, initial_guess, [], [], [], [], lb, ub, [], options);
% 提取结果
opt_a1 = best_params(1);
opt_a2 = best_params(2);
opt_b1 = best_params(3);
fprintf('拟合成功！\n最优 a1: %.4f\n最优 a2: %.4f\n最优 b1: %.4f\n', opt_a1, opt_a2, opt_b1);





%% --- 1. 结构基础与湿模态求解  ---
MM = m + Mw_global;
nm = 20;
[v,d] = eigs(k, MM, nm, 'SM');
tempd = diag(d);
[d, sortindex] = sort(tempd);
v = v(:, sortindex);
omega = sqrt(d);
frequency2 = omega / (2*pi);       % 转换频率(Hz)

f_ex = frequency2(1);              % 正确获取一阶湿模态频率作为激励频率
omega_ex = omega(1);               % 获取对应的圆频率
fprintf('当前激振频率 (一阶湿模态): %.4f Hz\n', f_ex);

%% --- 2. 阻尼矩阵组装 (必须在求出 omega_ex 之后) ---
zeta_1 = 0.010588;                 % 自由衰减实验的一阶阻尼比
beta_rayleigh = 2 * zeta_1 / omega_ex;
C_matrix = beta_rayleigh * k;      % 真实的结构粘性阻尼矩阵

% 组装总阻尼矩阵 (前提: Cw_global 已正确计算)
CC = C_matrix + Cw_global;

% 检查总阻尼矩阵是否因为 b1=-10.8 出现严重的负对角线元素
if min(diag(CC)) < 0
    warning('总阻尼矩阵存在负对角元素！流体呈现负阻尼，Newmark-Beta 必将发散！');
    warning('请重新运行 patternsearch，并务必设置 lb (下界) 强制 b1 > 0 !');
end

%% --- 3. 压电非线性饱和与迟滞等效参数 ---
alpha_sat = 1;                     % 300V下的极化饱和衰减系数
psi_hys_deg = 15;                  % 迟滞环导致的力学滞后角
psi_hys = psi_hys_deg * (pi / 180);% 转为弧度
h_act = 0.18e-3;                   % 实际真实的压电层厚度 (mm)
% 假设 hp 和 ka 已在前方定义
gamma_geo = h_act / hp ;
ka_true = ka * gamma_geo;
ka_eff = ka_true * alpha_sat;      % 修正机电耦合矩阵

%% --- 4. 时间历程与非线性激振力生成区 ---
V0 = -300;                       
Va = [V0; -V0];                    % 交流电压幅值向量 
Va_D = [0; 0];                     % 直流偏置电压向量 
phi = [0; 0];                      % 相位角差向量 

% 设定时间步和总时长
cycles = 24;                      % 仿真 100 个周期以确保达到稳态
T_end = cycles * (1 / f_ex);
points_per_cycle = 40;             % 每个周期采 40 个点保证精度
dt = (1 / f_ex) / points_per_cycle;
t = 0:dt:T_end;                    % 时间向量

% 在交流电压中引入物理滞后角 psi_hys，生成真实激振力矩阵
V_ac_eff = Va .* sin(phi + 2 * pi * f_ex * t - psi_hys); 
V_total_eff = V_ac_eff + Va_D;
F_matrix = ka_eff * V_total_eff;  

%% --- 5. 时域积分求解区 ---
fprintf('开始 Newmark-Beta 时域积分，共计算 %d 个时间步...\n', length(t));
[t_out1, X_out] = Standard_Newmark_Beta2(MM, CC, k, F_matrix, t);
fprintf('积分完成！\n');

% 绘制特定自由度的时域响应 
jth=8; % 观察的自由度8
time_vector = 0:dt:T_end;
figure(3)                   
plot(time_vector,X_out(jth,:)*1e3)
xlabel('Time (seconds)')
ylabel('displacement(mm)')
hold on
% ylim([-2.5, 2.5]);
xlim([min(time_vector), max(time_vector)]);

time_vector= 2.0:dt:2.7;
figure(4)                   
plot(time_vector,X_out(jth,2/dt:2.7/dt)*1e3);%401:801
xlabel('Time (seconds)')
ylabel('displacement(mm)')
hold on
xlim([min(time_vector), max(time_vector)]);


MM = m + Mw_global;           % 总质量
CC = C_matrix + Cw_global;    % 总阻尼
sys_dof = size(k, 1);              % 系统自由度数




% 1. 构建状态空间矩阵 A (维度 2N x 2N)
% 原方程: M*x'' + C*x' + K*x = 0
% 状态空间: X' = A*X, 其中 A = [0, I; -M^-1*K, -M^-1*C]
fprintf('正在组装状态空间矩阵...\n');
I = eye(sys_dof);
O = zeros(sys_dof);
A = [O, I; 
    -MM\k, -MM\CC];

% 2. 求解复特征值
fprintf('正在求解复特征值...\n');
% 注意：由于矩阵变为非对称，eig会计算出共轭复数根
[V_complex, D_complex] = eig(full(A)); 
lambda = diag(D_complex); % 提取特征值

% 3. 解析特征值
% 复特征值的形式为: lambda = -sigma +/- i * omega_d
% 其中虚部 omega_d 就是带阻尼的实际振荡圆频率
omega_d_all = abs(imag(lambda));

% 过滤掉静止解（虚部接近0），并且因为是共轭对，取唯一值
tol_freq = 1e-3;
omega_d_valid = unique(round(omega_d_all(omega_d_all > tol_freq), 4));

% 排序并转换为 Hz
freq_complex_Hz = sort(omega_d_valid) / (2*pi);

% 4. 打印对比结果
% 为了对比，顺手算一下你之前的不带阻尼的频率 (实模态)
[~, d_real] = eigs(k, MM, 1, 'SM');
freq_real_Hz = sqrt(d_real) / (2*pi);

fprintf('--------------------------------------------------\n');
fprintf('不带阻尼的实模态频率 (omega_n): %.4f Hz\n', freq_real_Hz);
fprintf('带总阻尼的复模态频率 (omega_d): %.4f Hz\n', freq_complex_Hz(1));
fprintf('两者频率差值: %.4f Hz\n', abs(freq_real_Hz - freq_complex_Hz(1)));
fprintf('--------------------------------------------------\n');