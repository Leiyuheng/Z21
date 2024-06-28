close all
clear
clc

c=3e8;
f = 3e9; %频率
k = 2*pi*f/c;
Rf = [1,0,1]*1e-3;
Rs = [0,0,1]*1e-3;
polar = [0,0,1]; 

%%============================首先计算线磁流产生的磁场====================
len = 2; % mm
dis = 20; % 中心点间距 mm
point = 64; % 采样点数
xs = linspace(-len/2,len/2,point)*1e-3;
xf = linspace(dis-len/2,dis+len/2,point)*1e-3;
E = zeros(point,point);
Js = sin(k*xs);
Jf = sin(k*(xf-dis*1e-3));

for i = 1:length(xs)
    for j = 1:length(xf)
        rs = [xs(i) Rs(2) Rs(3)];
        rf = [xf(j) Rf(2) Rf(3)];
        disp(rf)
        Ge = half_space_gf(rf,rs,f); % 场源相对位置（层数）确定 GF确定
        polar = [polar(1),polar(2),Js(i)];
        E(i,j) = calculate_E(Ge,polar,f);
    end
end
% E的每行代表每一个xs在不同xf下的场，每列代表不同xs在同一xf下的场
%%===============================积分计算==================
% 首先要将E的每一列对xs进行积分 压缩为一维
E_f = zeros(1,length(xf));
for i = 1:length(xf)
    E_s_re = real(E(i,:));
    E_s_im = imag(E(i,:));
    f_E_re = @(x) interp1(xs,E_s_re,x,'makima');
    f_E_im = @(x) interp1(xs,E_s_im,x,'makima');
    f_E = @(x) f_E_re(x) + 1i * f_E_im(x);

    E_f(i) = quadgk(f_E,xs(1),xs(end));
end

% 通过反应原理进行互阻抗计算
I1 = 1;
I2 = 1;
f_inte_disc = Jf.* E_f;
f_inte_re = @(x) interp1(xf,real(f_inte_disc),x,'makima');
f_inte_im = @(x) interp1(xf,imag(f_inte_disc),x,'makima');
f_inte = @(x) f_inte_re(x) + 1i * f_inte_im(x);

Z21 = quadgk(f_inte,xf(1),xf(end))/(I1*I2);
abs(Z21)
% % 打开文件
% fileID = fopen('90.efe', 'r');
% 
% % 读取文件头部信息，假设文件头部信息有固定行数
% headerLines = 15; % 跳过行数
% for i = 1:headerLines
%     fgetl(fileID); % 逐行读取并忽略文件头部信息
% end

% % 读取数据
% data = textscan(fileID, '%f %f %f %f %f %f %f %f %f', 'HeaderLines', 1);
% % 关闭文件
% fclose(fileID);
% 
% % 提取数据列
% x_feko = data{1};
% y_feko = data{2};
% z_feko = data{3};
% 
% Ex_feko = data{4} + 1i*data{5};
% Ey_feko = data{6} + 1i*data{7};
% Ez_feko = data{8} + 1i*data{9};
% E_feko = sqrt(Ex_feko.^2 + Ey_feko.^2 +Ez_feko.^2);
% 
% 
% % 绘制图形
% figure;
% hold on;
% 
% % 绘制仿真数据 (实线，绝对值)
% plot(z_feko, abs(Ex_feko), 'r-', 'DisplayName', 'Ex (仿真)');
% plot(z_feko, abs(Ey_feko), 'g-', 'DisplayName', 'Ey (仿真)');
% plot(z_feko, abs(Ez_feko), 'b-', 'DisplayName', 'Ez (仿真)');
% 
% % 绘制计算数据 (点，绝对值)
% plot(z_f, abs(Ex), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Ex (计算)');
% plot(z_f, abs(Ey), 'g*', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Ey (计算)');
% 
% % 设置图例和标签
% legend show;
% xlabel('Observation point x (m)');
% ylabel('Magnetic field (A/m)');
% title('Comparison of Simulated and Calculated Magnetic Fields');
% grid on;
% hold off;
