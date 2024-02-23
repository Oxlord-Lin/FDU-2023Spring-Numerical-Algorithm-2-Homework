t = [1,3,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]; %最后25使其成为完整的周期
y = 36 + [37,23,21,26,38,49,60,63,66,68,69,73,74,78,82,84,87,86,77,59,37]*1e-2; %最后一个也是为了插值方便而补充
x = [1:3:25];
h1 = scatter(t,y,'cyan','filled','s');
hold on
% only for test
% t = x;
% y = sin(t);
% test end

n = length(x) - 1;
A = zeros(3*n,4*n);
A(1:3,1:4) = funa(x(1));
A(1:3,end-3:end) = (-1)*funa(x(n+1)); % A是约束矩阵
for i = 2:n
    A((i-1)*3+1:(i)*3,(i-2)*4+1:(i-1)*4) = funa(x(i));
    A((i-1)*3+1:(i)*3,(i-1)*4+1:(i)*4) = (-1)*funa(x(i));
end
cla = zeros(8,1);
for i = 1:length(y) - 1
    cla(category(t(i))) = cla(category(t(i))) + 1; %得到落入每个区间的点的数量
end
cla(end) = cla(end) + 1; %最后一个点比较特殊

T = zeros(sum(cla),4*n);
T(1:cla(1),1:4) = funt(t(1:cla(1)));
for i = 2:8 %构造矩阵T
    T(sum(cla(1:i-1))+1:sum(cla(1:i)), 4*(i-1)+1:4*i ) = funt(t(sum(cla(1:i-1))+1:sum(cla(1:i))));
end
% 调个凸优化的包来处理，我实在不想写程序了
cvx_begin quiet
    variable g(32)
    r = y' - T*g; %r是残差
    minimize(norm(r))
    subject to 
        A*g == zeros(3*n,1);
cvx_end

% 绘图部分
for i = 1:n
    g_i = g((i-1)*4+1:i*4);
    a = g_i(1);
    b = g_i(2);
    c = g_i(3);
    d = g_i(4);
    tt = linspace(x(i),x(i+1),200);
    yy = zeros(200,1);
    f = @(u) a*u^3 + b*u^2 + c*u + d;
    for p = 1:200
        yy(p) = f(tt(p));
    end
    h = plot(tt,yy,tt+24,yy,Color='r',LineWidth=1);
    hold on
end
h2 = plot(tt,yy,Color='r',LineWidth=1); % 之后绘图方便
%% 
x = [1,3,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]; %最后25使其成为完整的周期
y = 36 + [37,23,21,26,38,49,60,63,66,68,69,73,74,78,82,84,87,86,77,59,37]*1e-2; %最后一个也是为了插值方便而补充
n = length(y);

% 准备工作
delta_x = zeros(n-1,1);
delta_y = zeros(n-1,1);
for i = 1:n-1
    delta_x(i) = x(i+1) - x(i);
    delta_y(i) = y(i+1) - y(i);
end
delta_x(n) = delta_x(1);
delta_x(n+1) = delta_x(2);
delta_y(n) = delta_y(1);
%
A = zeros(n);
A(1,1) = 1;
A(1,n) = -1;
for i = 2:n-1
    A(i,i) = 2*(delta_x(i)+delta_x(i-1));
    A(i,i-1) = delta_x(i);
    A(i,i+1) = delta_x(i-1);
end
A(n,1) = delta_x(n);
A(n,n-1) = delta_x(n+1);
A(n,n) = 2*(delta_x(n) + delta_x(n+1));
%
b = zeros(n,1);
for i = 2:n
    b(i) = 3*delta_x(i-1)*delta_y(i)/delta_x(i) + 3*delta_x(i)*delta_y(i-1)/delta_x(i-1);
end
k = A\b;
% 下面是绘图工作
for j = 1:n-1
    t = linspace(x(j),x(j+1),200);
    c = (delta_y(j)/delta_x(j) - k(j)) / delta_x(j); % Hermite插值的二阶系数
    d = (k(j+1) - 2*delta_y(j)/delta_x(j) + k(j))/(delta_x(j)^2); % Hermite插值的三阶系数
    f = zeros(200,1);
    for p = 1:200
        f(p) = y(j) + k(j)*(t(p)-x(j)) + c*((t(p)-x(j))^2) + d*((t(p)-x(j))^2)*(t(p)-x(j+1));
    end
    h = plot(t,f,t+24,f,LineWidth=1,LineStyle='--',Color='b');
    hold on
end
h3 = plot(t(1),f(1),LineWidth=1.2,LineStyle='--',Color='b'); % 之后绘图方便
ylim([36,37])
xlabel('time')
ylabel('averaged values of temperature')
title('Comparison of fitting and interpolating periodic cubic splines')

%%

legend([h1,h2,h3],'数据点','拟合','插值','Location','best')


%%
function A = funa(x)
    A = [x^3,x^2,x,1;
        3*x^2,2*x,1,0;
        6*x,2,0,0];
end

function T = funt(t)
    ki = length(t);
    T = ones(ki,4);
    for i = 1:ki
        T(i,1) = t(i)^3;
        T(i,2) = t(i)^2;
        T(i,3) = t(i);
    end
end

function cla = category(index)
    cla = floor((index-1)/3) + 1;
end