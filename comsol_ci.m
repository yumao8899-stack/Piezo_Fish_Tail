clc
clear
No_INTpoint_x=4;
No_INTpoint_y=4;
No_INTpoint_z=2;

% Bottom plate 层
% 线弹性材料每行为[E1, E2, G23, G12, G13, nu12, h, theta]
densityb=1620;

bottom_layers = [ 
    126e9, 8.73e9, 9.5e9, 7.8e9, 9.5e9, 0.31, 0.15e-3, 90];

middle_layers = [
    126e9, 8.73e9, 9.5e9, 7.8e9, 9.5e9, 0.31, 0.15e-3, 0;
    126e9, 8.73e9, 9.5e9, 7.8e9, 9.5e9, 0.31, 0.15e-3, 0;];

% Top plate 层
top_layers = [ 
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
hb=0.15e-3;
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

bp=[bt_p,bb_p];
ap=[at_p,ab_p];
ac=5.497;bc=3.65E-06;
% ac=14.206;bc=3.10E-06;
% ac=23.782;bc=3.52E-06;
%% 单元编号
[nodes, elements, elements_index] = read_comsol('yuwei.mphtxt');
elements = elements(:, [1 2 4 3]);
[nodes_ext, quads8] = convert_to_8node_shell(nodes, elements);
nodes_ext=1e-3*nodes_ext;
elements_p =quads8(elements_index{1},:);
%% 节点的物理坐标
% 上层节点坐标
node_total = size(nodes_ext,1);  %总节点数
coop = zeros(node_total, 3);  %初始化坐标矩阵，避免索引越界
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

% 节点平均坐标
% cooo = (coop + cop) / 2;
% x1 = coop(:,1); y1 = coop(:,2); z1 = coop(:,3);
% x2 = cop(:,1); y2 = cop(:,2); z2 = cop(:,3);

% 绘制初始节点位置
% figure(1)
% plot3(x1, y1, z1, 'b.'); hold on
% plot3(x2, y2, z2, 'r.'); hold off
% title('初始节点位置'); xlabel('X'); ylabel('Y'); zlabel('Z');

%%
%压电片的数量与位置
np=2;nn=elements_index{1};
sys_dof = dof;
k = zeros(sys_dof, sys_dof);    %系统刚度矩阵 
m = zeros(sys_dof, sys_dof);    %系统质量矩阵 
CC = zeros(sys_dof, sys_dof);


ka = zeros(sys_dof, np); 
ks = zeros(np, sys_dof); 
F = zeros(sys_dof, 1);
F(903,3) = -100;  
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
    ecpt = ac * empt + bc * ekpt;
    ecpb = ac * empb + bc * ekpb;
    ecpp = ac * empp + bc * ekpp;
    BB = ismember(loopi, nn);
    if BB ~= 0
        [ekp1,D] = middleshellek_p(ab_p,bb_p, middle_pian, H, dyhm, jdzbp, jdzb1p, dybh);
        [emm1] = middleshellem(ab_p,bb_p, zmtemp, v3i, densityp, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        [ekp2] = middleshellek_p(at_p,bt_p, middle_pian, H, dyhm, jdzbp, jdzb1p, dybh);
        [emm2] = middleshellem(at_p,bt_p, zmtemp, v3i, densityp, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        ecm1=ac * emm1 + bc * ekp1;
        ecm2=ac * emm2 + bc * ekp2;
        ek = ekpt + ekp1 + ekp2 + ekpb + ekpp;
        em = empt + emm1 + emm2 + empb + empp;
        ec = ecpt + ecm1 + ecm2 + ecpb + ecpp;
    else
        [eke1] = middleshellek(ab_p, bb_p, E_e, nu_e, H, dyhm, jdzbp, jdzb1p, dybh);
        [emm1] = middleshellem(ab_p, bb_p, zmtemp, v3i, densitye, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        [eke2] = middleshellek(at_p, bt_p, E_e, nu_e, H, dyhm, jdzbp, jdzb1p, dybh);
        [emm2] = middleshellem(at_p, bt_p, zmtemp, v3i, densitye, H, xv2i, xv1i, No_INTpoint_x, No_INTpoint_y, No_INTpoint_z);
        ecm1=ac * emm1 + bc * eke1;
        ecm2=ac * emm2 + bc * eke2;
        ek = ekpt + ekpb + eke2 + eke1 + ekpp;%
        em = empt + empb + emm2 + emm1 + empp;%
        ec = ecpt + ecm1 + ecm2 + ecpb + ecpp;
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
                CC(index(jx), index(jy)) = CC(index(jx), index(jy)) + ec(jx, jy);
            end
        end
    end
end

%% 组装压电耦合矩阵ka与ks
nnd=[nn,nn];
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

%% 模态分析
nm=20;
[v,d]=eigs(k,m,nm,'SM');
tempd = diag(d);
[d, sortindex] = sort(tempd);
v = v(:, sortindex);
omega=(sqrt(d));
frequency = (sqrt(d))/(2*pi);  %转换为频率(Hz)
% V1=[v(:,1:nm)];                                     
% Factor=diag(V1'*m*V1);
%% 模态形状可视化
node_total = size(disp_node,1);
nm = size(v,2);
v_full = zeros(node_total*5, nm);
map=disp_node';
map = map(:); % 拉平成一列
valid = map > 0;
v_full(valid, :) = v(map(valid), :);

mn=1;  %第1阶模态
voo = nodes_ext;
xvec = v_full;

% 提取各方向位移（假设自由度顺序为x,y,z）
xx = xvec(1:5:end, mn); 
yy = xvec(2:5:end, mn); 
zz = xvec(3:5:end, mn);  

voo(1:node_total , 1) = voo(1:node_total , 1) + 0.005*xx;
voo(1:node_total , 2) = voo(1:node_total , 2) + 0.005*yy;
voo(1:node_total , 3) = voo(1:node_total , 3) + 0.005*zz;

figure(2);  
title(['第 ', num2str(mn), ' 阶模态形状']);
faces4 = quads8(:, [1 2 3 4]);
patch('Faces', faces4, 'Vertices', voo, ...
      'FaceVertexCData',zz, ...
      'FaceColor', 'interp', 'EdgeColor', 'k');
hold off; 
%% 模态叠加法 
% 定义激励 
excitation_freq_Hz = 145.87; %[Hz] 
omega0 = excitation_freq_Hz * 2 * pi; % 转换为角频率
V0 = -600; 
Va =ones(np,1);
Va(1,1) = V0;Va(2,1) = -V0;
F0 =ka * Va; % 压电致动力 F = ka * V
F0(abs(F0)<1e-5)=0;
X_static = k \ F0;
[DIS,zz]= ttq(1,disp_node, X_static, nodes_ext,1);
figure(3);  
title(['第 ', num2str(mn), ' 阶模态形状']);
faces4 = quads8(:, [1 2 3 4]);
patch('Faces', faces4, 'Vertices', DIS, ...
      'FaceVertexCData',zz, ...
      'FaceColor', 'interp', 'EdgeColor', 'k');
hold off; 

dt = 5e-4; % 时间步长
tf = 2;  % 总时间
time_vector = 0:dt:tf;
C=eye(dof);q1=zeros(dof,1);dq0=zeros(dof,1);  
[eta,dip]=HarmonicRespt(k,m,F0,omega0,time_vector,C,q1,dq0,CC,omega,v);
% [eta,y,omega1,sdof2]=HarmonicRespt(k,m,F0,omega0,time_vector,C,q1,dq0,ac,bc);
jth=103;
figure(100)                   
plot(time_vector,dip(jth,:))
xlabel('Time (seconds)')
ylabel('displacement(m)')
xlim([min(time_vector), max(time_vector)]);

time_vector= 0.6:dt:0.8;
figure(102)                   
plot(time_vector,dip(jth,1201:1601))
xlabel('Time (seconds)')
ylabel('displacement(m)')
xlim([min(time_vector), max(time_vector)]);
ylim([-1.5e-3, 1.5e-3]);
hold on
%% 频响函数  
[f,P1]=FFT(time_vector,dip(jth,:));
figure(101)
plot(f, 20*log10(P1), 'b-','LineWidth', 1.5);
ylabel('Amplitude (dB)');
title('Frequency Response' );
xlabel('Frequency (Hz)');
grid on;
xlim([0, 1000]);
%% 变换到模态坐标
omega_n = frequency * 2 * pi; % 固有角频率向量
% 模态质量归一化 (v' * m * v = I)
for i = 1:nm
    v(:,i) = v(:,i) / sqrt(v(:,i)' * m * v(:,i));
end
% 计算模态阻尼比 (假设为瑞利阻尼)
zeta_n = (ac ./ (2*omega_n)) + (bc .* omega_n ./ 2);

% 计算模态力
f_n = v' * F0;

%% 求解模态的瞬态响应 q_n(t)
q_response = zeros(nm, length(time_vector)); % 存储所有模态的响应

for n = 1:nm
    % 每个模态的单自由度方程: q'' + 2ζω q' + ω² q = f*cos(ω₀t)
    omega = omega_n(n);
    zeta = zeta_n(n);
    force = f_n(n);
    
    % 使用解析解求解稳态+瞬态响应
    omega_d = omega * sqrt(1 - zeta^2); % 阻尼固有频率
    
    % 求解稳态部分
    H = 1 / (omega^2 - omega0^2 + 2i*zeta*omega*omega0);
    q_steady_state = force * abs(H);
    phi = angle(H);
    
    % 求解瞬态部分 (齐次解)，系数由初始条件 q(0)=0, q'(0)=0 确定
    A = -q_steady_state * cos(phi);
    B = (zeta*omega/omega_d) * A - (q_steady_state*omega0/omega_d) * sin(phi);
    
    % 叠加得到完整解
    for t_idx = 1:length(time_vector)
        t = time_vector(t_idx);
        q_transient = exp(-zeta*omega*t) * (A * cos(omega_d*t) + B * sin(omega_d*t));
        q_response(n, t_idx) = q_steady_state * cos(omega0*t + phi) + q_transient;
    end
end

%% 可视化
displacement_transient = v * q_response; % (dof x time_steps)
figure(99);
plot_dof = 103;
plot(time_vector, displacement_transient(plot_dof, :));
title('瞬态响应');
xlabel('时间 (s)');
ylabel('位移 (m)');
xlim([min(time_vector), max(time_vector)]);
%% 谐响应分析 
C = alpha_damping * m + beta_damping * k;
% 定义激励 
excitation_freq_Hz = 145.87; %[Hz] 
omega0 = excitation_freq_Hz * 2 * pi; % 转换为角频率

V0 = -600; % 假设施加的电压幅值为 100V
Va =ones(np,1);
Va(1,1) = V0;Va(2,1) = -V0;
F0 =ka * Va; % 压电致动力 F = ka * V

% 先使用原代码中的经验值
alpha_damping = 5; 
beta_damping = 0.000001;
C = alpha_damping * m + beta_damping * k;

D_matrix = k - omega0^2 * m + 1i * omega0 * C;
% 求解复数位移响应 X 
X_complex = D_matrix \ F0;
% 从复数响应中分离出同相(in-phase)和异相(out-of-phase)分量
% x(t) = real(X) * cos(ωt) + imag(X) * sin(ωt)
% 原代码中的 A0, B0 分别对应 imag(X) 和 real(X)

Xc = real(X_complex); % 与激励力同相的分量
Xs = -imag(X_complex);% 与激励力异相(滞后90度)的分量

% 计算位移幅值和相位 
amplitude = abs(X_complex);      % 每个自由度的振动幅值
phase = angle(X_complex);        % 每个自由度的相位角 (相对于激励力)
%%  后处理与可视化 
dt = 1 / (excitation_freq_Hz * 20); % 时间步长，保证每个周期至少有20个点
tf = 5 / excitation_freq_Hz;        % 总仿真时间，显示5个周期

tf = 0.25;
time_vector = 0:dt:tf;              % 创建时间向量

% 计算每个时间点的位移
displacement_over_time = zeros(dof, length(time_vector));
for i = 1:length(time_vector)
    t = time_vector(i);
    displacement_over_time(:, i) = Xc * cos(omega0 * t) + Xs * sin(omega0 * t);
end

% 绘制某个节点的位移-时间曲线 
figure(14);
plot_dof = 103; % 选择一个自由度来绘制
plot(time_vector, displacement_over_time(plot_dof, :)*1e3, 'b-', 'LineWidth', 1.5);
title(sprintf('自由度 %d 的位移-时间响应', plot_dof));
xlabel('时间 (s)');
ylabel('位移 (mm)');
grid on;

% 绘制位移最大的时刻的变形云图 
[~, max_disp_time_index] = max(abs(displacement_over_time(plot_dof, :)));
displacement_snapshot = displacement_over_time(:, max_disp_time_index);

[DIS,zz]= ttd(disp_node, displacement_snapshot, nodes_ext);
figure(15)
xlabel('X'); ylabel('Y'); zlabel('位移 (mm)');
faces4 = quads8(:, [1 2 3 4]);
patch('Faces', faces4, 'Vertices', DIS, ...
      'FaceVertexCData',zz, ...
      'FaceColor', 'interp', 'EdgeColor', 'k');
%title('位移最大的时刻变形');


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
