function [E_total,Ex,Ey,Ez] = calculate_E(Ge,polar,f)

%% ============================计算电场============================
c = 3e8;
u0 = 4*pi*1e-7; % 磁导率
e0 = 8.8542e-12; % 电导率
w = f * 2 * pi; % 角频率
% 首先定义电场 幅度为1 理想点源 极化方向为polar
% 由于要计算E的三个分量 所以将球坐标转为笛卡尔坐标


[Jx,Jy,Jz] = sph2cart(polar(1),polar(2),polar(3)); 
J = [Jx,Jy,Jz];
E = zeros(3,1);
E(:,:) = -1i * w * u0 * Ge * J';

Ex = squeeze(E(1, 1,:));
Ey = squeeze(E(2, 1,:));
Ez = squeeze(E(3, 1,:));
E_total = sqrt(Ex.^2 + Ey.^2 + Ez.^2);
end