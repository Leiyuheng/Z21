function [H_total,Hx,Hy,Hz] = calculate_H(Gh,polar,f)

%% ============================计算电场============================
c = 3e8;
u0 = 4*pi*1e-7; % 磁导率
e0 = 8.8542e-12; % 电导率
w = f * 2 * pi; % 角频率
% 首先定义电场 幅度为1 理想点源 极化方向为polar
% 由于要计算E的三个分量 所以将球坐标转为笛卡尔坐标


[Jx,Jy,Jz] = sph2cart(polar(1),polar(2),polar(3)); 
J = [Jx,Jy,Jz];
H = zeros(3,1);
H(:,:) = -1i * w * e0 * Gh * J';

Hx = squeeze(H(1, 1,:));
Hy = squeeze(H(2, 1,:));
Hz = squeeze(H(3, 1,:));
H_total = sqrt(Hx.^2 + Hy.^2 + Hz.^2);
end