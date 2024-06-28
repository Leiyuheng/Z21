function Gh = half_space_gf_cal(r, r_prime, f)
    % 计算半自由空间的并矢格林函数
    % 输入：
    %   R - 3x1 向量，观察点坐标 [x; y; z]
    %   R_prime - 3x1 向量，源点坐标 [x_prime; y_prime; z_prime]
    %   k - 波数
    % 输出：
    %   Ge1_val - 半自由空间的并矢格林函数值
    c = 3e8;
    lambda = c/f;
    k = 2*pi/lambda;
    % 定义符号变量
    syms x y z x_prime y_prime z_prime
    % 定义距离向量
    R = sqrt((x - x_prime)^2 + (y - y_prime)^2 + (z - z_prime)^2);
    Ri = sqrt((x - x_prime)^2 + (y - y_prime)^2 + (z + z_prime)^2);
    
    % 定义A, Ai, B, Bi
    A = exp(-1i*k*R)/(4*pi*R) * (3/(k^2 * R^2) + 3*1i/(k*R) - 1);
    Ai = exp(-1i*k*Ri)/(4*pi*Ri) * (3/(k^2 * Ri^2) + 3*1i/(k*Ri) - 1);
    B = exp(-1i*k*R)/(4*pi*R) * (1 - 1i/(k*R) - 1/(k^2 * R^2));
    Bi = exp(-1i*k*Ri)/(4*pi*Ri) * (1 - 1i/(k*Ri) - 1/(k^2 * Ri^2));
    
    % 定义格林函数分量
    Gxx = A/(R^2) * (x - x_prime)^2 + B + Ai/(Ri^2) * (x - x_prime)^2 + Bi;
    Gxy = A/(R^2) * (x - x_prime) * (y - y_prime) + Ai/(Ri^2) * (x - x_prime) * (y - y_prime);
    Gxz = A/(R^2) * (x - x_prime) * (z - z_prime) + Ai/(Ri^2) * (x - x_prime) * (z + z_prime);
    Gyy = A/(R^2) * (y - y_prime)^2 + B + Ai/(Ri^2) * (y - y_prime)^2 + Bi;
    Gyz = A/(R^2) * (y - y_prime) * (z - z_prime) + Ai/(Ri^2) * (y - y_prime) * (z + z_prime);
    Gzz = A/(R^2) * (z - z_prime)^2 + B - Ai/(Ri^2) * (z + z_prime)^2 - Bi;

%     curl_x = curl([Gxx, Gxy, Gxz], [x, y, z]);
%     curl_y = curl([Gxy, Gyy, Gyz], [x, y, z]);
%     curl_z = curl([Gxz, Gyz, Gzz], [x, y, z]);
% 
%     Gh = [curl_x';curl_y';curl_z'];
% 
    Gh = [Gxx, Gxy, Gxz; 
      Gxy, Gyy, Gyz; 
      Gxz, Gyz, Gzz];

    Gh = subs(Gh, [x, y, z, x_prime, y_prime, z_prime], [r(1), r(2), r(3), r_prime(1), r_prime(2), r_prime(3)]);
end
