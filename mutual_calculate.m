function Z21 = mutual_calculate(f,Rs_start,Rs_end,Rf_start,Rf_end,polar_s,polar_f,points)
% 参数意义介绍
% f: 工作频率
% Rs_start: 线磁流（源）的起始位置；Rs_end: 线磁流（源）的结束位置
% Rf_start: 线磁流（场）的起始位置；Rf_end: 线磁流（场）的结束位置
% polar_s：线磁流（源）的极化方向；polar_f：线磁流（场）的极化方向 应该是两个二维向量，不包含幅度，分别为phi和theta,强度由J决定
% points：采样点数

% Js: 源分布；Jf: 场分布 初始定义都是sin(kr)，其中k是波数，r是距离起始点距离； 如果要修改需要在本文件中修改

c = 3e8;
lambda = c/f;
k = 2*pi/lambda;
w = 2*pi*f;
e0 = 8.8542e-12; % 电导率
u0 = 4*pi*1e-7; % 磁导率
% 定义符号变量，r是场点，rs是源点
syms t % t是描述曲线的参数
rs = Rs_start + t * (Rs_end - Rs_start);

rf = [linspace(Rf_start(1),Rf_end(1),points)' linspace(Rf_start(2),Rf_end(2),points)' linspace(Rf_start(3),Rf_end(3),points)'];
Hf = zeros(points,3);

for i = 1:points
    R = norm(rf(i)-rs);
    g_R = exp(-1i*k*R)/(4*pi*R);
    A = g_R *(3/(k^2 * R^2)+3*1i/(k*R)-1);
    B = g_R *((1-1i/(k*R)-1/(k^2 * R^2)));
    
    rs_i = [rs(1),rs(2),-rs(3)];
    Ri = norm(rf(i)-rs_i);
    g_Ri = exp(-1i*k*Ri)/(4*pi*Ri);
    Ai = g_Ri *(3/(k^2 * Ri^2)+3*1i/(k*Ri)-1);
    Bi = g_Ri *((1-1i/(k*Ri)-1/(k^2 * Ri^2)));
    
    
    xf=rf(i,1);yf=rf(i,2);zf=rf(i,3);
    xs=rs(1);ys=rs(2);zs=rs(3);
    Js = sin(k*norm(rs-Rs_start));  % 源分布定义
    
    Ghxx = A/R^2 *(xf-xs)^2 + B + Ai/Ri^2 *(xf-xs)^2 + Bi;
    Ghxy = A/R^2 *(yf-ys)*(xf-xs) + Ai/Ri^2 *(yf-ys)*(xf-xs);
    Ghxz = A/R^2 *(zf-zs)*(xf-xs) + Ai/Ri^2 *(zf+zs)*(xf-xs);
    
    Ghyy = A/R^2 *(yf-ys)^2 + B + Ai/Ri^2 *(yf-ys)^2 + Bi;
    Ghyz = A/R^2 *(yf-ys)*(zf-zs) + Ai/Ri^2 *(yf-ys)*(zf+zs);
    
    Ghzz = A/R^2 *(zf-zs)^2 + B - Ai/Ri^2 *(zf+zs)^2 - Bi;
    
    Gh = [Ghxx Ghxy Ghxz; Ghxy Ghyy Ghyz; Ghxz Ghyz Ghzz];

    [Jx,Jy,Jz] = sph2cart(polar_s(1),polar_s(2),Js); 
    J = [Jx,Jy,Jz];
    H = -1i * w * e0 * Gh * J';

    % 将符号表达式转换为函数句柄
    H_func = matlabFunction(H, 'Vars', t);

    % 使用积分函数对 E 进行数值积分
    H_integral = integral(H_func, 0, 1, 'ArrayValued', true);
    Hf(i, :) = H_integral;
end

% 通过反应原理进行互阻抗计算
I1 = 1;
I2 = 1;
Hf_inte_disc = zeros(1,points);
for i = 1:points
    Jf = sin(k*norm(rf(i,:) - Rf_start)); % 场分布定义
    [Jx,Jy,Jz] = sph2cart(polar_f(1), polar_f(2),Jf);
    J = [Jx,Jy,Jz]';
    Hf_inte_disc(i) = Hf(i,:)*J;
end

r_int = zeros(1,points);
for i = 1:points
    r_int(i) = norm(rf(i,:) - Rf_start);
end
f_inte_re = @(x) interp1(r_int,real(Hf_inte_disc),x,'makima');
f_inte_im = @(x) interp1(r_int,imag(Hf_inte_disc),x,'makima');
f_inte = @(x) f_inte_re(x) + 1i * f_inte_im(x);

Z21 = quadgk(f_inte,r_int(1),r_int(end))/(I1*I2);


end




