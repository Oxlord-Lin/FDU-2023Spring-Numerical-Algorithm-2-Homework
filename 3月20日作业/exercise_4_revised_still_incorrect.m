%% 周期样条
x = [1,3,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,25]; %最后25使其成为完整的周期
y = 36 + [37,23,21,26,38,49,60,63,66,68,69,73,74,78,82,84,87,86,77,59,37]*1e-2; %最后一个也是为了插值方便而补充

% 测试部分
x = linspace(0,2*pi,4);
y = sin(x);
% 测试结束

n = length(y);
% 准备工作
delta_x = zeros(n-1,1);
delta_y = zeros(n-1,1);
% for i = 1:n-1
%     delta_x(i) = x(i+1) - x(i);
%     delta_y(i) = y(i+1) - y(i);
% end
for i = 1:n-1
    delta_x(i) = x(i+1) - x(i);
    delta_y(i) = y(i+1) - y(i);
end
%
A = zeros(n-1,1);
for j = 2:n-2
    i = j+1;
    A(j,j) = 2;
    A(j,j-1) = delta_x(i)/(delta_x(i)+delta_x(i-1));
    A(j,j+1) = delta_x(i-1)/(delta_x(i)+delta_x(i-1));
end
A(1,1) = 2;
A(1,2) = delta_x(1)/(delta_x(2) + delta_x(1));
A(1,n-1) = delta_x(2)/(delta_x(2) + delta_x(1));
A(n-1,1) = delta_x(n-1)/(delta_x(1) + delta_x(n-1));
A(n-1,n-2) = delta_x(1)/(delta_x(1) + delta_x(n-1));
A(n-1,n-1) = 2;
%
b = zeros(n-1,1);
for i = 1:n-2
    b(i) = 3 * (delta_x(i)/(delta_x(i)+delta_x(i+1))*(delta_y(i+1) / delta_x(i+1)) ...
        + delta_x(i+1)/(delta_x(i)+delta_x(i+1))*(delta_y(i) / delta_x(i)));
end
b(n-1) = 3 * (delta_x(n-1)/(delta_x(n-1)+delta_x(1))*(delta_y(1)/delta_x(1)) ...
    + delta_x(1)/(delta_x(n-1)+delta_x(1))*(delta_y(n-1)/delta_x(n-1)));

temp = A\b;
k = zeros(n-1,1);
k(1) = temp(n-1);
k(2:n) = temp(1:n-1);
% 下面是绘图工作
% figure
hold on
for j = 1:n-1
    t = linspace(x(j),x(j+1),200);
    c = (delta_y(j)/delta_x(j) - k(j)) / delta_x(j); % Hermite插值的二阶系数
    d = (k(j+1) - 2*delta_y(j)/delta_x(j) + k(j))/delta_x(j); % Hermite插值的三阶系数
    f = zeros(200,1);
    for p = 1:200
        f(p) = y(j) + k(j)*(t(p)-x(j)) + c*(t(p)-x(j))^2 + d*(t(p)-x(j))^2*(t(p)-x(j+1));
    end
    plot(t,f,LineWidth=1.5,Color='b');
    hold on
    plot(t+max(x),f,LineWidth=1.5,Color='c');
    hold on

end
plot(x,y,'o',Color='r')
xlabel('time')
ylabel('averaged values of temperature')
title('a periodic cubic spline of the two-day-period for the averaged temperature')