fun = @(x) x^64 - 0.1; % 原函数

%% regula falsi方法求根
a = 0;
y_a = fun(a);
b = 1;
y_b = fun(b);
c = (b*y_a - a*y_b)/(y_a - y_b);
regula_falsi_abs_er = zeros(1,1); % 用于记录收敛过程
regula_falsi_abs_er(1,1) = abs(fun(c));
count = 1;
while abs(fun(c))> 1e-12 % 使用abs(fun(c))作为误差判断依据
    count = count + 1;
    f_c = fun(c);
    if f_c == 0 
        break
    elseif f_c > 0
        b = c;
        f_b = f_c;
    else
        a = c;
        f_a = f_c;
    end
    c = (b*y_a - a*y_b)/(y_a - y_b);
    regula_falsi_abs_er(count,1) = abs(fun(c)); % 记录误差
end
semilogy(regula_falsi_abs_er,':',LineWidth=2)
hold on

%% 二分法求根
a = 0;
b = 1;
c = 0.5*(a+b);
bisection_abs_er = zeros(1,1);
bisection_abs_er(1,1) = abs(fun(c));
count2 = 1;
while abs(fun(c)) > 1e-12
    count2 = count2 + 1;
    if fun(c)>0
        b = c;
        c = (a+b)/2;
    elseif fun(c) == 0
        break
    else 
        a = c;
        c = (a+b)/2;
    end
    bisection_abs_er(count2,1) = abs(fun(c));
end
semilogy(bisection_abs_er,'--',LineWidth=2)
xlabel("迭代次数")
ylabel('误差（取对数）')
legend('absolute error by regula falsi','absolute error by bisection')
title('Bisection与regula falsi两种方法求方程 x^{64}-0.1=0 根的收敛过程')