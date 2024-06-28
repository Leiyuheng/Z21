% 基于三维格林函数与磁流的互耦计算

close all
clear
clc

c=3e8;
f = 3e9; %频率
k = 2*pi*f/c;

Rs_start = [-20,0,1.6]*1e-3;
Rs_end = [-30,0,1.6]*1e-3;
Rf_start = [-10,0,1.6]*1e-3;
Rf_end = [10,0,1.6]*1e-3;
polar_s = [0,0];
polar_f = [0,0]; % \theta and \phi
points = 64;

Z21 = zeros(1,points);
for i = 1:points+1
    phi = pi/points*(i-1);
    [x,y] = pol2cart(phi,10);
    Rf_start = [-x,-y,1.6]*1e-3;
    Rf_end = [x,y,1.6]*1e-3;
    polar_f = [phi,0];
    Z21(i) = mutual_calculate(f,Rs_start,Rs_end,Rf_start,Rf_end,polar_s,polar_f,points);
    
end

% 计算绝对值、实部和虚部
absZ = abs(Z21); % 绝对值
realZ = real(Z21); % 实部
imagZ = imag(Z21); % 虚部

% 创建时间向量或频率向量，假设与数据点数相同
phi = linspace(0, pi, length(Z21)); % 您可能需要根据实际情况调整这个向量

% 绘制数据
figure;
plot(phi, absZ, 'b-', 'LineWidth', 2); hold on; % 绝对值，蓝线
plot(phi, realZ, 'r--', 'LineWidth', 2); % 实部，红色虚线
plot(phi, imagZ, 'g-.', 'LineWidth', 2); % 虚部，绿色点划线
hold off;

% 添加图例
legend('Magnitude |Z|', 'Real Part Re(Z)', 'Imaginary Part Im(Z)', 'Location', 'best');

% 添加标题和坐标轴标签
title('Impedance Analysis');
xlabel('theta');
ylabel('Impedance Values');

xticks([0, pi/4, pi/2, 3*pi/4, pi]);
xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});

% 优化图形显示
grid on; % 打开网格



