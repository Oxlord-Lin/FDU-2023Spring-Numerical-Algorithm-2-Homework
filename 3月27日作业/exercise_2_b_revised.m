data = [3, 0; 2, 2; 0, 3; -2, 2; -3, 0; -2,-2; 0,-3;2,-2; 3,0]; % 最后自己再补一个点（即回到第一个点）
x = data(:,1);
y = data(:,2);
n = length(x);
t = zeros(n,1);
for i = 2:n
    t(i) = norm(data(i,:)-data(i-1,:)) + t(i-1);
end
% 自己写的周期样条插值有问题，之后需要再改（仍然不太对，不是饱满的圆）
xt = periodic_interpolate(t,x);
yt = periodic_interpolate(t,y);

% 用matlab内置的样条函数
% tt = 0:.0001:t(end);
% xt = spline(t,x,tt);
% yt = spline(t,y,tt);

plot(xt,yt,':',LineWidth=1)
axis equal
% hold on
% tt = 0:0.01:max(t);
% xx = spline(t,x,tt);
% yy = spline(t,y,tt);
% plot(xx,yy) % 这里用的是'not-a-knot'样条
hold on
plot(x,y,'*')



function f = periodic_interpolate(x,y)
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
count = 0;
f = zeros(200,1);
for j = 1:n-1
    t = linspace(x(j),x(j+1),200);
    c = (delta_y(j)/delta_x(j) - k(j)) / delta_x(j); % Hermite插值的二阶系数
    d = (k(j+1) - 2*delta_y(j)/delta_x(j) + k(j))/(delta_x(j)^2); % Hermite插值的三阶系数
    for p = 1:200
        count = count + 1;
        f(count) = y(j) + k(j)*(t(p)-x(j)) + c*((t(p)-x(j))^2) + d*((t(p)-x(j))^2)*(t(p)-x(j+1));
    end
end
end