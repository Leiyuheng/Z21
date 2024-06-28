function Ge = half_space_gf(r, r_prime, f)
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
    R = norm(r-r_prime);
    R_vec = (r-r_prime)./R;
    R_dyadic = [R_vec(1).*R_vec;R_vec(2).*R_vec;R_vec(3).*R_vec];
    I = eye(3);

    g_R = exp(-1i*k*R)/(4*pi*R);
    Ge0 = g_R *((3/(k^2 * R^2)+3*1i/(k*R)-1)*R_dyadic + (1-1i/(k*R)-1/(k^2 * R^2))*I);

    r_i_prime = [r_prime(1),r_prime(2),-r_prime(3)];
    R = norm(r-r_i_prime);
    R_vec = (r-r_i_prime)./R;
    R_dyadic = [R_vec(1).*R_vec;R_vec(2).*R_vec;R_vec(3).*R_vec];
    I = eye(3);

    g_R = exp(-1i*k*R)/(4*pi*R);
    Ge1 = g_R *((3/(k^2 * R^2)+3*1i/(k*R)-1)*R_dyadic + (1-1i/(k*R)-1/(k^2 * R^2))*I);

    nn_dyadic = [0 0 0; 0 0 0; 0 0 1];

    Ge = Ge0 - Ge1 + Ge1 .* nn_dyadic;
end
