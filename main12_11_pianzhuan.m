clc;clear
No_INTpoint_x=4;
No_INTpoint_y=4;
No_INTpoint_z=2;

% Bottom plate 层[E1, E2, G23, G12, G13, nu12, h, theta]
densityb=1620;
bottom_layers = [ 
    126e9, 8.73e9, 9.5e9, 7.8e9, 9.5e9, 0.31, 0.15e-3, 90
    126e9, 8.73e9, 9.5e9, 7.8e9, 9.5e9, 0.31, 0.15e-3, 45];

middle_layers = [
    126e9, 8.73e9, 9.5e9, 7.8e9, 9.5e9, 0.31, 0.15e-3, -45
    126e9, 8.73e9, 9.5e9, 7.8e9, 9.5e9, 0.31, 0.15e-3, -45];

top_layers = [ 
     126e9, 8.73e9, 9.5e9, 7.8e9, 9.5e9, 0.31, 0.15e-3, 45
     126e9, 8.73e9, 9.5e9, 7.8e9, 9.5e9, 0.31, 0.15e-3, 90];

% Epoxy and pizplate（各向同性）
densitye=1100;
E_e = 3e9; nu_e = 0.33; H_e = 0.3e-3;
G_e=E_e/(2*(1+nu_e));

densityp=7500;
% 压电材料性材料每行为[E1, E2, G, nu12, nu21，h]
middle_pian = [ 
        30.336e9, 15.857e9, 5.515e9, 0.31, 0.162 ,0.3e-3];
e31=-6.5;e32=-6.5;e24=0;e15=0;E33=1.3e-8;
 
hp=0.3e-3;
hb=0.3e-3;
hmid=0.3e-3;
H=2*hp+2*hb+hmid;
% layers_data_p=[bottom_layers;middle_layers;top_layers];
% layers_data_e=[bottom_layers;E_e,E_e,G_e,G_e,G_e,nu_e,H_e,0;top_layers];
% z0_piezo=get_neutral(layers_data_p);
% z0_epoxy=get_neutral(layers_data_e);
bt=1;                        at=(0.5*hmid+hp)/(0.5*H);
bt_p=(0.5*hmid+hp)/(0.5*H);  at_p=(0.5*hmid)/(0.5*H);
bm=(0.5*hmid)/(0.5*H);       am=-(0.5*hmid)/(0.5*H);
bb_p=-(0.5*hmid)/(0.5*H);    ab_p=-(0.5*hmid+hp)/(0.5*H);
bb=-(0.5*hmid+hp)/(0.5*H);   ab=-1;

bp=[bt_p,bb_p,bt_p,bb_p];
ap=[at_p,ab_p,at_p,ab_p];
%% 单元编号
[nodes, elements,elements_index] = read_comsol8('11.26wangge.bdf');
nodes_ext=1e-3*nodes;
combined_index = [elements_index{1},elements_index{3}];
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

%绘制初始节点位置
% figure(1)
% plot3(x1, y1, z1, 'b.'); hold on
% plot3(x2, y2, z2, 'r.'); hold off
% title('初始节点位置'); xlabel('X'); ylabel('Y'); zlabel('Z');
%% 压电片的数量与位置
np=4;
nn=combined_index;
sys_dof = dof;
k = zeros(sys_dof, sys_dof);    %系统刚度矩阵 
m = zeros(sys_dof, sys_dof);    %系统质量矩阵 
ka = zeros(sys_dof, np); 
ks = zeros(np, sys_dof); 
%%
jdzbp = coop;
jdzb1p = cop;
dybh = elements;
index = zeros(40, 1);  %存储单元节点的系统自由度
% 组装刚度矩阵和质量矩阵
for loopi=1:size(elements,1)
    dyhm = loopi;
    [ekpt] = topshellek_new(at, bt, top_layers, H, dyhm, jdzbp, jdzb1p, dybh);
    [ekpb, ~, xv2i, xv1i, ~, zmtemp, v3i, ~, jtemp] = bottomshellek_new(ab, bb, bottom_layers, H, dyhm, jdzbp, jdzb1p, dybh);
    [ekpp] = ppshellek_new(am, bm, middle_layers, H, dyhm, jdzbp, jdzb1p, dybh);
    [empt] = topshellem(at, bt, zmtemp, v3i, densityb, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
    [empb] = bottomshellem(ab, bb, zmtemp, v3i, densityb, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
    [empp] = ppshellem(am, bm, zmtemp, v3i, densityb, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
    BB = ismember(loopi, nn);
    if BB ~= 0
        [ekp1,D] = middleshellek_p(ab_p,bb_p, middle_pian, H, dyhm, jdzbp, jdzb1p, dybh);
        [emm1] = middleshellem(ab_p,bb_p, zmtemp, v3i, densityp, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        [ekp2] = middleshellek_p(at_p,bt_p, middle_pian, H, dyhm, jdzbp, jdzb1p, dybh);
        [emm2] = middleshellem(at_p,bt_p, zmtemp, v3i, densityp, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        ek = ekpt + ekp1 + ekp2 + ekpb + ekpp;
        em = empt + emm1 + emm2 + empb + empp;
    else
        [eke1] = middleshellek(ab_p, bb_p, E_e, nu_e, H, dyhm, jdzbp, jdzb1p, dybh);
        [emm1] = middleshellem(ab_p, bb_p, zmtemp, v3i, densitye, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        [eke2] = middleshellek(at_p, bt_p, E_e, nu_e, H, dyhm, jdzbp, jdzb1p, dybh);
        [emm2] = middleshellem(at_p, bt_p, zmtemp, v3i, densitye, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        ek = ekpt + ekpb + eke2 + eke1 + ekpp;%
        em = empt + empb + emm2 + emm1 + empp;%
    end
    for zi=1:8
        index((zi-1)*5 + 1) = disp_node(elements(loopi, zi), 1);
        index((zi-1)*5 + 2) = disp_node(elements(loopi, zi), 2);
        index((zi-1)*5 + 3) = disp_node(elements(loopi, zi), 3);
        index((zi-1)*5 + 4) = disp_node(elements(loopi, zi), 4);
        index((zi-1)*5 + 5) = disp_node(elements(loopi, zi), 5);
    end
% 组装整体矩阵
    for jx=1:40
        for jy=1:40
            if index(jx) ~= 0 && index(jy) ~= 0
                k(index(jx), index(jy)) = k(index(jx), index(jy)) + ek(jx, jy);
                m(index(jx), index(jy)) = m(index(jx), index(jy)) + em(jx, jy);
            end
        end
    end
end
%% 组装压电耦合矩阵ka与ks
nnd=[nn(:,1),nn(:,1),nn(:,2),nn(:,2)];
for n=1:np
    for loopi=1:size(elements,1)
        EE = ismember(loopi, nnd(:,n));
        if EE ~= 0
            dyhm = loopi;
            [ekat] = PZT_topka(ap(1,n), bp(1,n), e24, e15, H, hp, e31, e32, dyhm, jdzbp, jdzb1p, dybh);
          % [ekat] = MFC_topka(ap(1,n), bp(1,n), e24, e15, H, hp, e31, e32, dyhm, jdzbp, jdzb1p, dybh);
            eka = 2*ekat;
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
frequency = (sqrt(d))/(2*pi);  %转换为频率(Hz)

%% 模态叠加法 
% 定义激励 (简谐激励)[左上，左下，右上，右下]
excitation_freq_Hz=300;
excitation_full = [300,300,300,300]; %[Hz]； 102.6  % 4个不同的激励频率%102.6,102.6,102.6,102.6
omega0 = excitation_full * 2 * pi; % 转换为角频率
V0 = -600; 
Va = [V0; -V0; V0; -V0];
F0 =ka .* Va'; 
phi = [0; 0; pi; pi];
F0(abs(F0)<1e-5)=0;                                       

%直流分量
Va_D = [0; 0;0;0];%-600; 0; -600; 0        0; -600; 0; -600
F1 =ka * Va_D; 
F_full=[F0,F1];omega0_full=[omega0,0];  phi_full = [phi; 0];
% X_static = k \ F1;
% [DIS_D,dip_D]= ttq(disp_node, X_static, nodes_ext);
% mn=1;
% figure(3);
% title(['第 ', num2str(mn), ' 阶模态形状']);
% faces4 = elements(:, [1 2 3 4]);
% patch('Faces', faces4, 'Vertices', DIS, ...
%       'FaceVertexCData',zz, ...
%       'FaceColor', 'interp', 'EdgeColor', 'k');
% hold off; 

% 定义阻尼矩阵 (Rayleigh Damping)
% 测定阻尼比计算得到阻尼系数 (2%)
% zeta1 = 0.02; zeta2 = 0.02;
% w1 = omega(1); w2 = omega(2);
% A = 0.5 * [1/w1, w1; 1/w2, w2];
% b_damp = [zeta1; zeta2];
% alpha_beta = A \ b_damp;
% ac = alpha_beta(1);  % m 矩阵的系数
% bc = alpha_beta(2);  % k 矩阵的系数

ac=55.5; bc=3.65E-06;   %5.497
CC = ac * m + bc * k;
dt = 1/(excitation_freq_Hz*16);              %时间步长5e-4
tf = 2;                 %总时间
time_vector = 0:dt:tf;
C=eye(dof);q1=zeros(dof,1);dq0=zeros(dof,1);  
[eta,dip]=anti_HarmonicRespt(k,m,F_full,omega0_full,time_vector,C,q1,dq0,CC,omega,v,phi_full);

% 绘制特定自由度的时域响应
jth=3; % 观察的自由度3
time_vector = 0:dt:tf;
figure(100)                   
plot(time_vector,dip(jth,:)*1e3)
xlabel('Time (seconds)')
ylabel('displacement(mm)')
hold on
% ylim([-2.5, 2.5]);
xlim([min(time_vector), max(time_vector)]);

time_vector= 0.2:dt:0.4;
figure(102)                   
plot(time_vector,dip(jth,0.2/dt:0.4/dt)*1e3)%401:801
xlabel('Time (seconds)')
ylabel('displacement(mm)')
hold on
xlim([min(time_vector), max(time_vector)]);

[~, max_disp] = max(abs(dip(jth, 801:3777)));
q_snapshot = dip(:, 800+max_disp)*1e3;              %这一时刻的完整自由度位移向量
[DIS, zz] = ttd8(disp_node, q_snapshot, nodes_ext);  %zz建议用竖向位移或位移模长
figure(15)
scale_factor = 150; max_abs_z = max(abs(zz)); 
DIS_scaled = DIS; DIS_scaled(:, 3) = zz * scale_factor; 
figure(15); 
set(gcf, 'Color', 'w'); % 背景设为白色
faces4 = elements(:, 1:4);
hPatch = patch('Faces', faces4, ...
               'Vertices', DIS_scaled, ...    % 用放大的坐标画图
               'FaceVertexCData', abs(zz), ...  % 用真实值上色
               'FaceColor', 'interp', ...
               'EdgeColor', 'none');          % 无网格线，显示平滑云图
view(-30,40); box on; grid on; axis tight; axis equal;
set(gca, 'GridLineStyle', '--', 'LineWidth', 0.8, 'FontName', 'Times New Roman');
ylim_current = ylim;
set(gca, 'YTick', -40:20:40); % 设置一个足够大的范围，MATLAB会自动截取
num_ticks = 3; 
real_tick_values = linspace(-max_abs_z, max_abs_z, num_ticks); 
scaled_tick_positions = real_tick_values * scale_factor;
set(gca, 'ZTick', scaled_tick_positions); 
set(gca, 'ZTickLabel', num2str(real_tick_values', '%.2f')); % 保留3位小数
z_limit_view = max_abs_z * scale_factor;
zlim([-z_limit_view, z_limit_view]);
colormap("turbo"); 
clim([0, max_abs_z]); % 颜色范围 0 到 最大值
cb = colorbar;
title(cb, '(mm)', 'FontSize', 11, 'FontWeight', 'normal');
cb.TickLabelInterpreter = 'tex';