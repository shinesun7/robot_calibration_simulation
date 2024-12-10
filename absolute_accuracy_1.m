clear all;
clc;
% 符号计算
syms theta d a alpha beta d_theta d_d d_a d_alpha d_beta
T = trotx(alpha)*transl(a,0,0)*trotz(theta)*transl(0,0,d);%相邻坐标系齐次转换矩阵，MDH模型，T=RTTR

syms theta1 d1 a1 alpha1 beta1 theta2 d2 a2 alpha2 beta2 theta3 d3 a3 alpha3 beta3  
syms theta4 d4 a4 alpha4 beta4 theta5 d5 a5 alpha5 beta5 theta6 d6 a6 alpha6 beta6

T01 = subs(T, {alpha, a, theta, d}, {alpha1, a1, theta1, d1});      
T12 = subs(T, {alpha, a, theta, d}, {alpha2, a2, theta2, d2});
T23 = subs(T, {alpha, a, theta, d}, {alpha3, a3, theta3, d3});
T34 = subs(T, {alpha, a, theta, d}, {alpha4, a4, theta4, d4});
T45 = subs(T, {alpha, a, theta, d}, {alpha5, a5, theta5, d5});
T56 = subs(T, {alpha, a, theta, d}, {alpha6, a6, theta6, d6});

T06 = T01*T12*T23*T34*T45*T56;

params = [alpha1 a1 theta1 d1 alpha2 a2 theta2 d2 alpha3 a3 theta3 d3,...
          alpha4 a4 theta4 d4 alpha5 a5 theta5 d5 alpha6 a6 theta6 d6,];

% 直接T06对4*6个参数求导,如果使用24个参数的雅可比矩阵，
% 使用这个J_P,需要将J_p = zeros(3 * size(data, 1)/2, 20)中的20改为24
J_P = jacobian(T06(1:3,4),params);

% 去了冗余参数列d3、d5、theta5、theta6
J_P(:,[12,19,20,23]) = [];%3行20列

% 数值计算
% 名义DH参数
alpha = [  0, 90,   0, -90, 90,-90]*pi/180;
    a = [  0, 166.605, -782.270, -138.826,  0,  0];
    d = [504,   0,   0, 761.350,  0, 125];
theta_init = [0 -90 0 0 0 0]*pi/180;

data = readmatrix('LF.xlsx');
theta = data*pi/180 + theta_init;
% 仿真-实际DH参数
theta_new = theta+[0.05 -0.14 0.04 -0.17 0.24 -0.09]*pi/180;
    d_new = d+    [-0.42 0.36 -0.71 0.96 -0.32 0.38];
    a_new = a+    [0.73 0.54 -0.41 0.30 0.84 -0.05];
alpha_new = alpha+[0.21 -0.13 0.08 -0.16 0.37 0.11]*pi/180;

% 计算仿真实际末端位置
end_position_r = zeros(3*size(data,1)/2,1);%3*50行1列
for i=1:size(data,1)/2
    T01_r=subs(trotx(alpha_new(1))*transl(a_new(1),0,0)*trotz(theta_new(i,1))*transl(0,0,d_new(1)));
    T12_r=subs(trotx(alpha_new(2))*transl(a_new(2),0,0)*trotz(theta_new(i,2))*transl(0,0,d_new(2)));
    T23_r=subs(trotx(alpha_new(3))*transl(a_new(3),0,0)*trotz(theta_new(i,3))*transl(0,0,d_new(3)));
    T34_r=subs(trotx(alpha_new(4))*transl(a_new(4),0,0)*trotz(theta_new(i,4))*transl(0,0,d_new(4)));
    T45_r=subs(trotx(alpha_new(5))*transl(a_new(5),0,0)*trotz(theta_new(i,5))*transl(0,0,d_new(5)));
    T56_r=subs(trotx(alpha_new(6))*transl(a_new(6),0,0)*trotz(theta_new(i,6))*transl(0,0,d_new(6)));
    
    T06_r = T01_r*T12_r*T23_r*T34_r*T45_r*T56_r*[0;0;0;1];
    end_position_r((i-1)*3+1:i*3, :) = T06_r(1:3);
end

lambda = 0.0004;  % 初始 lambda 值
Maxiter = 10;  % 最大迭代次数

for iter = 1:Maxiter
    disp(['Iteration ', num2str(iter)]);
    disp(['d(1)=',num2str(d(1))]);

    end_position_n = zeros(3 * size(data, 1) / 2, 1);  % 估计的末端位置
    J_p = zeros(3 * size(data, 1) / 2, 20);  % 位置误差雅可比矩阵
    % 计算末端位置
    for i = 1:size(data, 1) / 2
        % 计算各个变换矩阵
        T01_n = subs(trotx(alpha(1)) * transl(a(1), 0, 0) * trotz(theta(i, 1)) * transl(0, 0, d(1)));
        T12_n = subs(trotx(alpha(2)) * transl(a(2), 0, 0) * trotz(theta(i, 2)) * transl(0, 0, d(2)));
        T23_n = subs(trotx(alpha(3)) * transl(a(3), 0, 0) * trotz(theta(i, 3)) * transl(0, 0, d(3)));
        T34_n = subs(trotx(alpha(4)) * transl(a(4), 0, 0) * trotz(theta(i, 4)) * transl(0, 0, d(4)));
        T45_n = subs(trotx(alpha(5)) * transl(a(5), 0, 0) * trotz(theta(i, 5)) * transl(0, 0, d(5)));
        T56_n = subs(trotx(alpha(6)) * transl(a(6), 0, 0) * trotz(theta(i, 6)) * transl(0, 0, d(6)));

        % 计算末端位置
        T06_n = T01_n * T12_n * T23_n * T34_n * T45_n * T56_n * [0; 0; 0; 1];
        end_position_n((i - 1) * 3 + 1:i * 3, :) = T06_n(1:3);

        % 设置参数
        params_n = [alpha(1) a(1) theta(i, 1) d(1) alpha(2) a(2) theta(i, 2) d(2), ...
                    alpha(3) a(3) theta(i, 3) d(3) alpha(4) a(4) theta(i, 4) d(4), ...
                    alpha(5) a(5) theta(i, 5) d(5) alpha(6) a(6) theta(i, 6) d(6)];

        % 计算、更新雅可比矩阵
        J_p((i - 1) * 3 + 1:i * 3, :) = double(subs(J_P, params, params_n));
    end

    % 计算末端位置的误差
    delta_P = end_position_r - end_position_n;%delta_x,delta_y,delta_z

    % 计算 LM 更新增量
    JtJ = J_p' * J_p;
    JtP = J_p' * delta_P;
    H = JtJ + lambda * eye(size(JtJ));  % 添加正则化项
    delta_X = H \ JtP;  % 使用 LM 算法求解增量
    
    rs = J_p * delta_X - delta_P;
    Rs = sum(rs.^2);
    disp(['Rs=',num2str(Rs)]);

    if Rs > 0.001
%         % 更新参数,
%         %alpha1 a1 theta1 d1 alpha2 a2 theta2 d2 alpha3 a3 theta3 0,...
%         %alpha4 a4 theta4 d4 alpha5 a5 0 0 alpha6 a6 0 d6,
        alpha = alpha + [delta_X(1) delta_X(5) delta_X(9) delta_X(12) delta_X(16) delta_X(18)];
            a = a + [delta_X(2) delta_X(6) delta_X(10) delta_X(13) delta_X(17) delta_X(19)];
            d = d + [delta_X(4) delta_X(8) 0 delta_X(15) 0 delta_X(20)];
        theta = theta + [delta_X(3) delta_X(7) delta_X(11) delta_X(14) 0 0]; 
    else
        break;
    end

end

error_alpha = (alpha - [  0, 90,   0, -90, 90, -90]*pi/180)*180/pi;
    error_a = a - [  0, 166.605, -782.270, -138.826,  0,  0];
error_theta = (theta(1,:) - [0 -90 0 0 0 0]*pi/180)*180/pi - data(1,:);
    error_d = d - [504,   0,   0, 761.350,  0, 125];

