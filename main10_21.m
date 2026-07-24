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
[nodes, elements, elements_index] = read_comsol4('10.31yu.mphtxt');
elements = elements(:, [1 2 4 3]);
[nodes_ext, quads8] = convert_to_8node_shell(nodes, elements);
nodes_ext=1e-3*nodes_ext;
combined_index = [elements_index{3},elements_index{7};elements_index{4}, elements_index{8}]; 
elements_p = quads8(combined_index, :);
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
zuobiaox = find(abs(nodes_ext(:,1) - 0.11) < 1e-6);
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
%nn=elements_index{1};
nn=combined_index;
sys_dof = dof;
k = zeros(sys_dof, sys_dof);    %系统刚度矩阵 
m = zeros(sys_dof, sys_dof);    %系统质量矩阵 
ka = zeros(sys_dof, np); 
ks = zeros(np, sys_dof); 
%%
jdzbp = coop;
jdzb1p = cop;
dybh = quads8;
index = zeros(40, 1);  %存储单元节点的系统自由度
% 组装刚度矩阵和质量矩阵
for loopi=1:size(quads8,1)
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
        index((zi-1)*5 + 1) = disp_node(quads8(loopi, zi), 1);
        index((zi-1)*5 + 2) = disp_node(quads8(loopi, zi), 2);
        index((zi-1)*5 + 3) = disp_node(quads8(loopi, zi), 3);
        index((zi-1)*5 + 4) = disp_node(quads8(loopi, zi), 4);
        index((zi-1)*5 + 5) = disp_node(quads8(loopi, zi), 5);
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
    for loopi=1:size(quads8,1)
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
            index((zi-1)*5 + 1) = disp_node(quads8(loopi, zi), 1);
            index((zi-1)*5 + 2) = disp_node(quads8(loopi, zi), 2);
            index((zi-1)*5 + 3) = disp_node(quads8(loopi, zi), 3);
            index((zi-1)*5 + 4) = disp_node(quads8(loopi, zi), 4);
            index((zi-1)*5 + 5) = disp_node(quads8(loopi, zi), 5);
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
% V1=[v(:,1:nm)];                                     
% Factor=diag(V1'*m*V1);
% %% 模态分析(water)
% nm=20;
% [v,d]=eigs(k,m+,nm,'SM');
% tempd = diag(d);
% [d, sortindex] = sort(tempd);
% v = v(:, sortindex);
% omega=(sqrt(d));
% frequency1 = (sqrt(d))/(2*pi);  %转换为频率(Hz)
%% 模态形状可视化
% node_total = size(disp_node,1);
% nm = size(v,2);
% v_full = zeros(node_total*5, nm);
% map=disp_node';
% map = map(:); % 拉平成一列
% valid = map > 0;
% v_full(valid, :) = v(map(valid), :);
% 
% mn=3;  %第1阶模态
% voo = nodes_ext;
% xvec = v_full;
% 
% % 提取各方向位移（假设自由度顺序为x,y,z）
% xx = xvec(1:5:end, mn); 
% yy = xvec(2:5:end, mn); 
% zz = xvec(3:5:end, mn);  
% 
% voo(1:node_total , 1) = voo(1:node_total , 1) + 0.005*xx;
% voo(1:node_total , 2) = voo(1:node_total , 2) + 0.005*yy;
% voo(1:node_total , 3) = voo(1:node_total , 3) + 0.005*zz;
% 
% figure(2);  
% title(['第 ', num2str(mn), ' 阶模态形状']);
% faces4 = quads8(:, [1 2 3 4]);
% patch('Faces', faces4, 'Vertices', voo, ...
%       'FaceVertexCData',abs(zz), ...
%       'FaceColor', 'interp', 'EdgeColor', 'none');
% colormap(turbo)
% cl = max(abs(zz(:)));
% if cl>0, clim([0 cl]); end
% hold off;
% axis equal;   % 保持正确的长宽比(colorbar显示图例)
%axis off;    % 关闭坐标轴、刻度和背景框
%% 模态叠加法 
% 定义激励 (简谐激励)[左上，左下，右上，右下]
excitation_freq_Hz=300;
excitation_full = [300,300,300,300]; %[Hz]； 102.6  % 4个不同的激励频率%102.6,102.6,102.6,102.6
omega0 = excitation_full * 2 * pi; % 转换为角频率
V0 = -600; 
Va = [V0; -V0; V0; -V0];%2*V0; -V0; 2*V0; -V0     V0; -2*V0; V0; -2*V0
F0 =ka .* Va'; 
phi = [0; 0; pi; pi];
F0(abs(F0)<1e-5)=0;                                      

%直流分量
Va_D = [0; 0; 0; 0];%-600; 0; -600; 0        0; -600; 0; -600
F1 =ka * Va_D; 
% X_static = k \ F1;
% [DIS_D,dip_D]= ttq(disp_node, X_static, nodes_ext);
F_full=[F0,F1];omega0_full=[omega0,0];  phi_full = [phi; 0];


% figure(3);
% title(['第 ', num2str(mn), ' 阶模态形状']);
% faces4 = quads8(:, [1 2 3 4]);
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
jth=18; % 观察的自由度23
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

[~, max_disp] = max(abs(dip(jth, 801:4001)));
q_snapshot = dip(:, 800+max_disp)*1e3;              %这一时刻的完整自由度位移向量
[DIS, zz] = ttd(disp_node, q_snapshot, nodes_ext);  %zz建议用竖向位移或位移模长
figure(15)
faces4 = quads8(:, [1 2 3 4]);
patch('Faces', faces4, 'Vertices', DIS, ...
'FaceVertexCData', abs(zz), ...
'FaceColor', 'interp', 'EdgeColor', 'none');
axis equal; colorbar; colormap(turbo)
title('挠曲变形图')

%% 遗传算法优化 
% x = [前段压电片, 后段压电片, 相位差, 直流分量]
lb = [0,  0,  0,  0];  % 下界
ub = [1000, 1000, pi,  0];  % 上界

%% 2. 目标选择（三选一）
% 目标 A: 最大化推力（直线游动）
fitness_func = @(x) -objective_thrust(x, k, m, ka, omega, v, disp_node, nodes_ext);

% 目标 B: 实现左转（转向角最大）
% fitness_func = @(x) -objective_turn(x, ...);

% 目标 C: 匹配目标振型
% fitness_func = @(x) -objective_shape(x, ..., target_shape);

%% 3. GA 设置
options = optimoptions('ga', ...
    'PopulationSize', 30, ...
    'MaxGenerations', 50, ...
    'Display', 'iter', ...
    'PlotFcn', @gaplotbestf);

%% 4. 运行优化
[x_opt, fval] = ga(fitness_func, 4, [], [], [], [], lb, ub, [], options);

%% 5. 输出结果
fprintf('\n========== 优化结果 ==========\n');
fprintf('前段电压: %.0f V\n', x_opt(1));
fprintf('后段电压: %.0f V\n', x_opt(2));
fprintf('相位差: %.1f °\n', x_opt(3)*180/pi);
fprintf('DC偏置: %.0f V\n', x_opt(4));
fprintf('适应度: %.4f\n', -fval);



%% 频响函数  
%frf
[f, U_amp] = calculate_frf(k, m, ka, Va, jth, omega, v, ac, bc, nm);
figure(101);
plot(f, U_amp*1e3, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Displacement Amplitude (mm)');
title('Sweep Frequency Response');
xlim([0 1000]);

%FFT
% [f,P1]=FFT(time_vector,dip(jth,:));
% figure(101)
% plot(f, 20*log10(P1), 'b-','LineWidth', 1.5);
% ylabel('Amplitude (dB)');
% title('Frequency Response' );
% xlabel('Frequency (Hz)');
% grid on;
% xlim([0, 200]);

%% 偏转
% [XD,YD]= tu_1(disp_node , dip, nodes_ext,excitation_freq_Hz);
% figure(500); 
% set(gcf, 'Color', 'w'); % 设置图形背景为白色
% hold on; % 允许在同一张图上叠加绘制
%     colors = [
%         0.8500 0.3250 0.0980; % 橙色 (T/4) 0.3010 0.7450 0.9330; % 浅蓝
%         0.9290 0.6940 0.1250; % 黄色
%         0.4660 0.6740 0.1880; % 绿色
%         0.4940 0.1840 0.5560; % 紫色
%         0   0.4470    0.7410; % 蓝色 (T/2, 中间时刻)
%         0.4940 0.1840 0.5560; % 紫色 
%         0.4660 0.6740 0.1880; % 绿色
%         0.9290 0.6940 0.1250; % 黄色
%         0.8500 0.3250 0.0980  % 橙色 (3T/4) 0.3010 0.7450 0.9330; % 浅蓝
%     ];
%     markers = {'o', '*', 'd', '^', 's', '^', 'd', '*', 'o'};
%     % 创建图例标签 (这是一个示例，你可以根据你的周期自己定义)
%     snapshot_labels = arrayfun(@(x) sprintf('t_{%d}', x), 1:9, 'UniformOutput', false);
% for i = 1:9
%     % 如果标记不够，则循环使用
%     marker_idx = mod(i-1, length(markers)) + 1;
%     plot(XD, YD(:, i), ...
%          'LineStyle', '-', ...
%          'Marker', markers{marker_idx}, ...
%          'Color', colors(i, :), ...
%          'LineWidth', 1.5, ...
%          'MarkerSize', 6, ...
%          'DisplayName', snapshot_labels{i}); % DisplayName 用于图例
% end
% hold off;
%% 验证
figure(12); hold on; axis equal;
% 画四节点单元的网格线quads8
for i = 1:size(elements,1)
    idx = elements(i, :); % 四节点单元
    xy = nodes_ext(idx, 1:2);
    plot([xy([1 2 3 4 1],1)], [xy([1 2 3 4 1],2)], 'k-');
end
% 只标出单元用到的节点编号
used_nodes = unique(elements(:));
for k = 1:length(used_nodes)
    i = used_nodes(k);
    text(nodes_ext(i,1), nodes_ext(i,2), num2str(i), ...
        'FontSize',12,'Color','b','HorizontalAlignment','center');
end
title('四节点单元网格及用到的节点编号');
xlabel('X'); ylabel('Y');


figure(11); hold on; axis equal;
for i = 1:size(quads8,1)
    idx = quads8(i,:); % [n1 n2 n3 n4 m12 m23 m34 m41]
    xy = nodes_ext(idx, 1:2);
    % 真实八节点壳单元边
    edge_order = [1 5 2 6 3 7 4 8 1]; % n1-m12-n2-m23-n3-m34-n4-m41-n1
    plot(xy(edge_order,1), xy(edge_order,2), 'k-');
    for j = 1:8
        text(xy(j,1), xy(j,2), num2str(idx(j)), ...
            'FontSize',12,'Color','b','HorizontalAlignment','center');
    end
end

figure(10); hold on; axis equal;
for i = 1:size(quads8,1)
    idx = quads8(i, [1 2 3 4]); % 八节点单元的四个角点
    xy = nodes_ext(idx, 1:2);
    % 画四条边
    plot([xy([1 2],1); NaN], [xy([1 2],2); NaN], 'k-');
    plot([xy([2 3],1); NaN], [xy([2 3],2); NaN], 'k-');
    plot([xy([3 4],1); NaN], [xy([3 4],2); NaN], 'k-');
    plot([xy([4 1],1); NaN], [xy([4 1],2); NaN], 'k-');
end
axis off
title('所有八节点单元角点边界线');
xlabel('X'); ylabel('Y');
