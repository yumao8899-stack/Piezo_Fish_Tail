clc;
clear;
%matlab_file = load('C:\Users\Administrator\Desktop\8.30\matlab位移\位移9.2.mat'); 
matlab_file = load('C:\Users\Administrator\Desktop\8.30\matlab位移\1.mat'); 

comsol_file = load('C:\Users\Administrator\Desktop\8.30\comsol位移\共振0—2.txt');
% time = comsol_file(401:2001, 1);
% disp1 = comsol_file(401:2001, 2);

time = comsol_file(:, 1);
disp1 = comsol_file(:, 2);

dof= 38; 
disp2 = -matlab_file.dip(dof,:)'*-1e3;
dispo = disp2;%;(1197:1597,1)

% 计算 TRAC 
numerator = (disp1' * dispo)^2;
denominator = (disp1' * dispo) * (disp1' * dispo);
trac_value = numerator / denominator;
% 计算MAE
absolute_errors = abs(disp1 - dispo);
mae_value = mean(absolute_errors);

% 绘制对比图
figure;
plot(time, disp1, 'b-', 'LineWidth', 2,'DisplayName', 'COMSOL ' );%
hold on;
plot(time, dispo, 'r--', 'LineWidth', 1.5,'DisplayName', 'matlab');% 
hold off;

grid on;

title(sprintf('TRAC = %.4f, MAE = %.4f', trac_value, 1.1397), 'FontSize', 14);
legend('show', 'Location', 'northeast');



%xlim([0.6, 0.8])



% xlabel('时间 (s)', 'FontSize', 12);
% ylabel('位移 (mm)', 'FontSize', 12);