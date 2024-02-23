data = [223,404
    906,138
    1410,38
    1264,191
    948,310
    1596,222
    1766,325
    1141,471
    1768,454
    1900,586
    1309,667
    1151,706
    1609,759
    1889,842
    1373,951
    910,970
    1269,1117
    1121,1248
    478,1144
    257,1063
    87,554
    223,404
    ];
x = data(:,1);
y = -1*data(:,2);
n = length(x);
t = zeros(n,1);
for i = 2:n
    t(i) = norm(data(i,:)-data(i-1,:)) + t(i-1);
end
% 自己写的周期样条插值有点问题，之后需要再改
% xt = periodic_interpolate(t,x);
% yt = periodic_interpolate(t,y);
% plot(xt,yt,LineWidth=1)
% hold on
tt = 0:0.01:max(t);
xx = spline(t,x,tt);
yy = spline(t,y,tt);
plot(xx,yy)
legend
hold on
plot(x,y,'*')
title('the 2D interpolation of the nodes taken from the outlier of my left hand')

% 周期样条的代码我从上一次作业搬过来的 但是有问题！写的不对
function yy = periodic_interpolate(x,y)
n = length(y);
delta_x = zeros(n-1,1);
delta_y = zeros(n-1,1);
for i = 1:n-1
    delta_x(i) = x(i+1) - x(i);
    delta_y(i) = y(i+1) - y(i);
end
%
A = zeros(n-1,1);
for j = 2:n-2
    i = j+1;
    A(j,j) = 2*(delta_x(i)+delta_x(i-1));
    A(j,j-1) = delta_x(i);
    A(j,j+1) = delta_x(i-1);
end
A(1,1) = delta_x(2) + delta_x(1);
A(1,2) = delta_x(1);
A(1,n-1) = delta_x(2);
A(n-1,1) = delta_x(n-1);
A(n-1,n-2) = delta_x(1);
A(n-1,n-1) = delta_x(1) + delta_x(n-1);
%
b = zeros(n-1,1);
for i = 1:n-2
    b(i) = 3 * (delta_x(i)*(delta_y(i+1) / delta_x(i+1)) ...
        + delta_x(i+1)*(delta_y(i) / delta_x(i)));
end
b(n-1) = 3 * (delta_x(n-1)*(delta_y(1)/delta_x(1)) ...
    + delta_x(1)*(delta_y(n-1)/delta_x(n-1)));

temp = A\b;
k = zeros(n-1,1);
k(1) = temp(n-1);
k(2:n) = temp(1:n-1);
yy = zeros(n,1);
count = 0;
for j = 1:n-1
    t = linspace(x(j),x(j+1),200);
    c = (delta_y(j)/delta_x(j) - k(j)) / delta_x(j); % Hermite插值的二阶系数
    d = (k(j+1) - 2*delta_y(j)/delta_x(j) + k(j))/delta_x(j); % Hermite插值的三阶系数
    for p = 1:200
        count = count + 1;
        yy(count) = y(j) + k(j)*(t(p)-x(j)) + c*(t(p)-x(j))^2 + d*(t(p)-x(j))^2*(t(p)-x(j+1));
    end
end
end


