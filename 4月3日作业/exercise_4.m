data=[
-4.00000 0.00001
-3.50000 0.00726
-3.00000 0.25811
-2.50000 1.87629
-2.00000 1.55654
-1.50000 0.17209
-1.00000 0.00899
-0.50000 0.05511
0.00000 0.24564
0.50000 0.60455
1.00000 0.89370
1.50000 1.03315
2.00000 0.51633
2.50000 0.18032
3.00000 0.04287
3.50000 0.00360
4.00000 0.00045];
x = data(:,1);
y = data(:,2);

c = [2,1,2,1,-2.5,1.5]'; % 做出良好的初始估计需要一些观察与估计
res = residual(x,y,c(1),c(2),c(3),c(4),c(5),c(6)); % 初始残差
J = Jacobi(x,c(1),c(2),c(3),c(4),c(5),c(6));
step = 1; % 步长
while step > 1e-15  % 我不太理解的是，为什么不能用norm(J)<1e-16 作为收敛的判别标准？
    c_new = c - pinv(J)*res; % 伪逆
    step = norm(c_new-c);
    c = c_new;
    res = residual(x,y,c(1),c(2),c(3),c(4),c(5),c(6));
    J = Jacobi(x,c(1),c(2),c(3),c(4),c(5),c(6));
end
xx = linspace(-4,4,1000);
yy = zeros(1000,1);
a1 = c(1);
a2 = c(2);
b1 = c(3);
b2 = c(4);
r1 = c(5);
r2 = c(6);
for i = 1:1000
    yy(i) = a1*exp(-(b1*(xx(i)-r1))^2) + a2*exp(-(b2*(xx(i)-r2))^2);
end
plot(xx,yy);
hold on
plot(x,y,'*')
title(['\alpha_1 =',num2str(a1),', \alpha_2=',num2str(a2),' \beta_1=',num2str(b1),' \beta_2=',num2str(b2),' \gamma_1=',num2str(r1),' \gamma_2=',num2str(r2)])
legend('拟合曲线','采样点')
function f = fun(x,a1,a2,b1,b2,r1,r2) % 计算某点的函数值
    f = a1*exp(-(b1*(x-r1))^2) + a2*exp(-(b2*(x-r2))^2);
end

function res = residual(x,y,a1,a2,b1,b2,r1,r2) % 计算某点的残差
    n = length(x);
    res = zeros(n,1);
    for i = 1:n
        res(i) = fun(x(i),a1,a2,b1,b2,r1,r2) - y(i);
    end
end

function J = Jacobi(x,a1,a2,b1,b2,r1,r2) % 计算某点（六个参数是自变量）的Jacobi矩阵
    f1 = @(t) exp(-(b1*(t-r1))^2);
    f2 = @(t) exp(-(b2*(t-r2))^2);
    f3 = @(t) a1*exp(-(b1*(t-r1))^2)*(t-r1)^2*(-2*b1);
    f4 = @(t) a2*exp(-(b2*(t-r2))^2)*(t-r2)^2*(-2*b2);
    f5 = @(t) a1*exp(-(b1*(t-r1))^2)*(-b1^2)*(2*(r1-t));
    f6 = @(t) a2*exp(-(b2*(t-r2))^2)*(-b2^2)*(2*(r2-t));
    n = length(x);
    J = zeros(n,6);
    for i = 1:n
        J(i,1) = f1(x(i));
        J(i,2) = f2(x(i));
        J(i,3) = f3(x(i));
        J(i,4) = f4(x(i));
        J(i,5) = f5(x(i));
        J(i,6) = f6(x(i));
    end
end