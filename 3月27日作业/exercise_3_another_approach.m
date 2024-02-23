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
semilogy(regula_falsi_abs_er,':',LineWidth=1,Color='g')
xlabel("迭代次数")
ylabel('误差（取对数）')
hold on
y = log10(regula_falsi_abs_er); % 直接用原始数据会出现条件数很大的情况，因此先取对数，再还原回去
n = length(y);
A = [(1:n)', ones(n,1)];
t = A\y;
k = t(1); % 斜率
b = t(2); % 截距
xx = 0:0.01:n;
yy = k*xx + b;
for i = 1:length(yy)
    yy(i) = 10^(yy(i));
end
hold on
plot(xx,yy,LineWidth=1,Color='b')
legend('the history of residual using regula falsi','最小二乘法拟合的指数类曲线（图中取了对数）')
[~,~] = title(['the history of residual using regula falsi and the curve that fits the history -- another approach'], ...
    ['k=',num2str(k),' b=',num2str(b)]);




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
figure
semilogy(bisection_abs_er,'--',LineWidth=1,Color='g')
xlabel("迭代次数")
ylabel('误差（取对数）')
hold on
% 选取残差曲线上的下述诸点，用最小二乘法进行直线拟合
y = log10(bisection_abs_er); % 直接用原始数据会出现条件数很大的情况，因此先取对数，再还原回去
n = length(y);
A = [(1:n)', ones(n,1)];
t = A\y;
k = t(1); % 斜率
b = t(2); % 截距
xx = 0:0.01:n;
yy = k*xx + b;
for i = 1:length(yy)
    yy(i) = 10^(yy(i));
end
plot(xx,yy,LineWidth=1,Color='b')
legend('the history of residual using bisection', '最小二乘法拟合的指数类曲线（图中取了对数）')
[~,~] = title(['the history of residual using bisection and the curve that fits the history -- another approach'], ...
    ['k=',num2str(k),' b=',num2str(b)]);