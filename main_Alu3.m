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
% d31=-1.7e-10;d33=4e-10;     %(1000vpp左右)
% d31=-1.4e-10; d33=3e-10;    %(500vpp)
d31=-1.3e-10; d33=2.51e-10; %(200vpp)


% Epoxy and pizplate（各向同性）
densitye=1100;
E_e = 3e9; nu_e = 0.32;
G_e=E_e/(2*(1+nu_e));


%damping
zeta=0.01378;
eta_b = 2*zeta;  % 铝板的损耗因子
eta_p = 2*zeta;  % 压电材料的损耗因子 

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

%%
jdzbp = coop;
jdzb1p = cop;
dybh = elements;
index = zeros(40, 1);  %存储单元节点的系统自由度
% 组装刚度矩阵和质量矩阵

for loopi = 1:size(elements,1)
    isMFC = ismember(loopi, nn);
    [ek_b, ~, xv2i, xv1i, ~, zmtemp, v3i, ~, jtemp] =middleshellek(am, bm, E_b, nu_b, H, loopi, jdzbp, jdzb1p, dybh);
    [em_b] = middleshellem(am, bm, zmtemp, v3i, densityb, H, xv2i, xv1i,No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
    ek_total = ek_b;
    em_total = em_b;
    ek_complex_local = ek_b * (1 + 1i * eta_b);
    if isMFC
        [ek_m1, D] = middleshellek_p(at, bt, middle_pian, H, loopi,jdzbp, jdzb1p, dybh);
        [em_m1] = middleshellem(at, bt, zmtemp, v3i, densityp, H,xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        [ek_m2] = middleshellek_p(ab, bb, middle_pian, H, loopi,jdzbp, jdzb1p, dybh);
        [em_m2] = middleshellem(ab, bb, zmtemp, v3i, densityp, H,xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        [ek_e1, ~] = middleshellek(a1, b1, E_e, nu_e, H, loopi,jdzbp, jdzb1p, dybh);
        [em_e1] = middleshellem(a1, b1, zmtemp, v3i, densitye, H,xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        [ek_e2, ~] = middleshellek(a2, b2, E_e, nu_e, H, loopi,jdzbp, jdzb1p, dybh);
        [em_e2] = middleshellem(a2, b2, zmtemp, v3i, densitye, H,xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);


        ek_total = ek_total + ek_m1 + ek_m2 + ek_e1+ ek_e2;
        em_total = em_total + em_m1 + em_m2 + em_e1+ em_e2;

        % 阻尼（MFC区域）
        ek_complex_local = ek_total * (1 + 1i * eta_p);

    end
    for zi = 1:8
        index((zi-1)*5 + 1) = disp_node(elements(loopi, zi), 1);
        index((zi-1)*5 + 2) = disp_node(elements(loopi, zi), 2);
        index((zi-1)*5 + 3) = disp_node(elements(loopi, zi), 3);
        index((zi-1)*5 + 4) = disp_node(elements(loopi, zi), 4);
        index((zi-1)*5 + 5) = disp_node(elements(loopi, zi), 5);
    end

    nonzero_idx = index(index ~= 0);
    valid_mask = index ~= 0;

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



%% --- 实验标定参数输入区 ---
% 1. 仅使用刚度比例阻尼 (根据自由衰减实验的一阶衰减结果填入)
f_ex = frequency(1);  
omega_ex = 2 * pi * f_ex;              % 圆频率
beta_rayleigh = 2 * zeta / omega_ex;
C_matrix1 = beta_rayleigh * k;         % 真实的粘性阻尼矩阵

C_matrix2 = imag(K_complex) / omega_ex; %复模态阻尼



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
Va = [V0; 0];                    % 交流电压幅值向量 [压电片1; 压电片2]
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



F0=4*F_matrix2;
X_static =k \F0;
[DIS_D,dip_D]= ttq(disp_node, X_static, nodes_ext);
figure(11)
faces4 = elements(:, [1 2 3 4]);
patch('Faces', faces4, 'Vertices', DIS_D, ...
'FaceVertexCData', abs(dip_D), ...
'FaceColor', 'interp', 'EdgeColor', 'none');
axis equal; colorbar; colormap(turbo)
title('挠曲变形图')
%% 模态叠加法 
% 定义激励 (简谐激励)[左上，左下，右上，右下]
omega0_full=[omega_ex,0];
phi_full=[phi;0];

cycles = 100;                    % 仿真 100 个周期以确保达到稳态
T_end = cycles * (1 / f_ex);
points_per_cycle = 40;           % 每个周期采 40 个点保证精度
dt = (1 / f_ex) / points_per_cycle;
t = 0:dt:T_end;                  % 时间向量 (尺寸: 1 x Nt)                       
C=eye(dof);q1=zeros(dof,1);dq0=zeros(dof,1);  
[eta,dip]=anti_HarmonicRespt(k,m,F_matrix2,omega0_full,t,C,q1,dq0,C_matrix2,omega,v,phi);

% 绘制特定自由度的时域响应
jth=8; % 观察的自由度8
time_vector = 0:dt:T_end;
figure(1)                   
plot(time_vector,dip(jth,:)*1e3)
xlabel('Time (seconds)')
ylabel('displacement(mm)')
hold on
% ylim([-2.5, 2.5]);
xlim([min(time_vector), max(time_vector)]);

time_vector= 2.0:dt:2.6;
figure(2)                   
plot(time_vector,dip(jth,2.0/dt:2.6/dt)*1e3);
xlabel('Time (seconds)')
ylabel('displacement(mm)')
hold on
xlim([min(time_vector), max(time_vector)]);



%% 频域求解器 
jth=8;
f_vec = 1 : 0.5 : 100; 
[f_out, X_complex, X_amp, X_phase] = Harmonic_Response_Solver(m, K_complex, F_matrix2, f_vec);

figure(5);
% 绘制幅值曲线 (通常采用对数坐标更能看清共振峰)
subplot(2,1,1);
semilogy(f_out, X_amp(jth, :), 'b-', 'LineWidth', 1.5);
grid on; hold on;
title(sprintf('DOF %d 的幅频特性曲线', jth));
xlabel('频率 (Hz)');
ylabel('位移幅值 (m)');
xlim([min(f_out), max(f_out)]);

% 绘制相位曲线
subplot(2,1,2);
plot(f_out, X_phase(jth, :) * 180/pi, 'r-', 'LineWidth', 1.5);
grid on; hold on;
title(sprintf('DOF %d 的相频特性曲线', jth));
xlabel('频率 (Hz)');
ylabel('相位角 (Degree)');
xlim([min(f_out), max(f_out)]);
ylim([-180, 180]);
yticks(-180:90:180);
