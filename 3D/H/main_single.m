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
Rf = [0 1 1.5]*1e-3;
zf = linspace(-1,1,point)*1e-3;
H = zeros(1,point);
Hx = zeros(1,point);
Hy = zeros(1,point);
Hz = zeros(1,point);

for i = 1:point
    rf = [zf(i) Rf(2) Rf(3)];
    Gh = half_space_gf_cal(rf,rs,f); % 场源相对位置（层数）确定 GF确定
    [H(i),Hx(i),Hy(i),Hz(i)] = calculate_H(Gh,polar,f);
end


% 打开文件
fileID = fopen('90.hfe', 'r');

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

Hx_feko = data{4} + 1i*data{5};
Hy_feko = data{6} + 1i*data{7};
Hz_feko = data{8} + 1i*data{9};
H_feko = sqrt(Hx_feko.^2 + Hy_feko.^2 +Hz_feko.^2);


% 绘制图形
figure;
hold on;

% 绘制仿真数据 (实线，绝对值)
plot(x_feko, abs(H_feko), 'r-', 'DisplayName', 'Hx (仿真)');
% plot(x_feko, abs(Ey_feko), 'g-', 'DisplayName', 'Ey (仿真)');
% plot(x_feko, abs(Ez_feko), 'b-', 'DisplayName', 'Ez (仿真)');

% 绘制计算数据 (点，绝对值)
plot(zf, abs(H), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Hx (计算)');
% plot(xf, abs(Ey), 'g*', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Ey (计算)');
% plot(xf, abs(Ez), 'b*', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Ez (计算)');

% 设置图例和标签
legend show;
xlabel('Observation point x (m)');
ylabel('Magnetic field (A/m)');
title('Comparison of Simulated and Calculated Magnetic Fields');
grid on;
hold off;
