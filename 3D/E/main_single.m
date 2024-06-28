close all
clear
clc

c=3e8;
f = 3e9; %频率
k = 2*pi*f/c;

%%============================首先计算线磁流产生的磁场====================

point = 64; % 采样点数
polar = [0,0,1];
rs =[1 1 1]*1e-3;
Rf = [0 1 1]*1e-3;
zf = linspace(0.1,2.1,point)*1e-3;
E = zeros(1,point);
Ex = zeros(1,point);
Ey = zeros(1,point);
Ez = zeros(1,point);

for i = 1:point
    rf = [Rf(1) Rf(2) zf(i)];
    Ge = half_space_gf_cal(rf,rs,f); % 场源相对位置（层数）确定 GF确定
    [E(i),Ex(i),Ey(i),Ez(i)] = calculate_E(Ge,polar,f   );
end


% 打开文件
fileID = fopen('90.efe', 'r');

% 读取文件头部信息，假设文件头部信息有固定行数
headerLines = 15; % 跳过行数
for i = 1:headerLines
    fgetl(fileID); % 逐行读取并忽略文件头部信息
end
% 读取数据
data = textscan(fileID, '%f %f %f %f %f %f %f %f %f', 'HeaderLines', 1);
% 关闭文件
fclose(fileID);

% 提取数据列
x_feko = data{1};
y_feko = data{2};
z_feko = data{3};

Ex_feko = data{4} + 1i*data{5};
Ey_feko = data{6} + 1i*data{7};
Ez_feko = data{8} + 1i*data{9};
E_feko = sqrt(Ex_feko.^2 + Ey_feko.^2 +Ez_feko.^2);


% 绘制图形
figure;
hold on;

% 绘制仿真数据 (实线，绝对值)
plot(z_feko, abs(E_feko), 'r-', 'DisplayName', 'Ex (仿真)');
% plot(x_feko, abs(Ey_feko), 'g-', 'DisplayName', 'Ey (仿真)');
% plot(x_feko, abs(Ez_feko), 'b-', 'DisplayName', 'Ez (仿真)');

% 绘制计算数据 (点，绝对值)
plot(zf, abs(E), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Ex (计算)');
% plot(xf, abs(Ey), 'g*', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Ey (计算)');
% plot(xf, abs(Ez), 'b*', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Ez (计算)');

% 设置图例和标签
legend show;
xlabel('Observation point x (m)');
ylabel('Magnetic field (A/m)');
title('Comparison of Simulated and Calculated Magnetic Fields');
grid on;
hold off;
