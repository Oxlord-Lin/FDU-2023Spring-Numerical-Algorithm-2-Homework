%% 用"Good" & "Bad" Broyeden method求解雨水的pH值
%% 准备部分
% 各平衡常数都抽象成常量
c1 = 1e-14;
c2 = 10^(-6.3 - 6 - 1.46) * 375;
c3 = 10^(-10.3);
% 用于记录误差变化情况
er1 = zeros(1,1); % 记录Good Broyden的收敛过程
er2 = zeros(1,1); % 记录Bad Broyden的收敛过程
% 初始状态
x = 1e-7;  % x是氢离子浓度
y = 1e-7;  % y是氢氧根浓度
z = 1e-7;  % z是碳酸氢根浓度
w = 1e-7;  % w是碳酸根浓度
v_0 = [x,y,z,w]';  % 将各离子浓度存储为列向量v
F = myfunc(v_0);
er1(1,1) = norm(F); 
er2(1,1) = er1(1,1);
%% "Good" Broyden
% 得到的解一直不太正确
v = v_0; 
J = eye(4);
inv_J = eye(4);
u = zeros(4,1);
w = u;
F = myfunc(v);
count = 1;
while norm(F) > 1e-11 %精度无法继续提高了
    count = count + 1;

    % 以下利用Shermann-Morrison-Woodbury公式求逆
    inv_J_new = inv_J - inv_J*u*w'*inv_J./(1 + w'*inv_J*u); % 这里的inv_J是上一轮得到的
    inv_J = inv_J_new; 

    norm(inv_J*J-eye(4),'fro') % 值得注意的是，通过S-M-W公式得到的inv_J质量不好，却仍然可以收敛
    N(count-1,1) = norm(inv_J*J-eye(4),'fro');

    v2 = v - inv_J*myfunc(v); % 更新v

    v2 = v - J\myfunc(v); % 若用常规的方式更新v，但效果不一定会更好

    delta_v = v2 - v;
    F2 = myfunc(v2);
    delta_F = F2 - F;
    
    u = (delta_F - J*delta_v)./norm(delta_v);
    w = delta_v./norm(delta_v);

    J2 = J + u*w';
%      J2 = ((delta_F - J*delta_v)*delta_v')./(delta_v'*delta_v);

    v = v2;
    F = F2;
    J = J2;
 

    er1(count,1) = norm(F);
end
semilogy(er1)





%% 计算F
function F = myfunc(v)
c1 = 1e-14;
c2 = 10^(-6.3 - 6 - 1.46) * 375;
c3 = 10^(-10.3);
x = v(1,1);
y = v(2,1);
z = v(3,1);
w = v(4,1);
F = [x*y - c1;
    x*z - c2;
    x*w/z - c3;
    x-y-z-2*w;];
end




