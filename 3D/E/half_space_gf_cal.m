function Ge = half_space_gf_cal(r, r_prime, f)
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
    g_R = exp(-1i*k*R)/(4*pi*R);
    A = g_R *(3/(k^2 * R^2)+3*1i/(k*R)-1);
    B = g_R *((1-1i/(k*R)-1/(k^2 * R^2)));

    r_i_prime = [r_prime(1),r_prime(2),-r_prime(3)];
    Ri = norm(r-r_i_prime);
    g_Ri = exp(-1i*k*Ri)/(4*pi*Ri);
    Ai = g_Ri *(3/(k^2 * Ri^2)+3*1i/(k*Ri)-1);
    Bi = g_Ri *((1-1i/(k*Ri)-1/(k^2 * Ri^2)));


    xf=r(1);yf=r(2);zf=r(3);
    xs=r_prime(1);ys=r_prime(2);zs=r_prime(3);

    Gexx = A/R^2 *(xf-xs)^2 + B - Ai/Ri^2 *(xf-xs)^2 -Bi;
    Gexy = A/R^2 *(yf-ys)*(xf-xs) - Ai/Ri^2 *(yf-ys)*(xf-xs);
    Gexz = A/R^2 *(zf-zs)*(xf-xs) - Ai/Ri^2 *(zf+zs)*(xf-xs);

    Geyy = A/R^2 *(yf-ys)^2 + B - Ai/Ri^2 *(yf-ys)^2 -Bi;
    Geyz = A/R^2 *(yf-ys)*(zf-zs) - Ai/Ri^2 *(yf-ys)*(zf+zs);

    Gezz = A/R^2 *(zf-zs)^2 + B + Ai/Ri^2 *(zf+zs)^2 +Bi;

    Ge = [Gexx Gexy Gexz; Gexy Geyy Geyz; Gexz Geyz Gezz];
end
