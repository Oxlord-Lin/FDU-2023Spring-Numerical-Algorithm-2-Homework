n = 9;
x = linspace(0,2*pi,n);
y = sin(x);
delta_x = zeros(8);
delta_y = delta_x;
for i = 1:n-1
    delta_x(i) = x(i+1) - x(i);
    delta_y(i) = y(i+1) - y(i);
end

A = diag(ones(n,1));
for i = 2:n-1
    A(i,i) = 2*(delta_x(i)+delta_x(i-1));
    A(i,i-1) = delta_x(i);
    A(i,i+1) = delta_x(i-1);
end
b = zeros(n,1);
for i = 2:n-1
    b(i) = 3*delta_x(i-1)*(delta_y(i)/delta_x(i)) + 3*delta_x(i)*(delta_y(i-1)/delta_x(i-1));
end
b(1) = NaN; % 待定
b(n) = 1;


slope = [-4,-2,4,2]; % x = 0 处的一阶导
for i = 1:length(slope)
    b(1) = slope(i); %一阶导
    k = A\b;
    % 下面是绘图工作
    ax = nexttile;
    for j = 1:n-1
        t = linspace(x(j),x(j+1),200);
        c = (delta_y(j)/delta_x(j) - k(j))/delta_x(j); % Hermite插值的二阶系数
        d = (k(j+1) - 2*(delta_y(j)/delta_x(j)) + k(j))/delta_x(j)^2; % Hermite插值的三阶系数
        f = zeros(200,1);
        for p = 1:200
            f(p) = y(j) + k(j)*(t(p)-x(j)) + c*(t(p)-x(j))^2 + d*(t(p)-x(j))^2*(t(p)-x(j+1));
        end
        Colors = ['g','c','r','m'];
        col = Colors(i); % 每一段的颜色都要统一
        plot(ax,t,f,Color=col,LineWidth=1.5);
        hold on
    end
    temp_x = linspace(0,2*pi,1000);
    temp_y = sin(temp_x);
    plot(ax,temp_x,temp_y,'k--',LineWidth=1);
    title(['x=0处的一阶导为',num2str(slope(i)),'，虚线为sin(x)图像'])
    ylim([-1.1,1.1])
    xlim([0,2*pi])
end