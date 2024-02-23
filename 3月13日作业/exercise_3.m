nodes = [2,3,5,9];
for i = 1:4
    ax = nexttile;
    t = linspace(0,2*pi,1000);
    f = sin(t);
    plot(t,f,'--',LineWidth=1.5,Color='r');
    hold on
    node_num = nodes(i);
    x = linspace(0,2*pi,node_num);
    y = sin(x);
    y_diff = cos(x);
    cubic_Hermite_interpolation(x,y,y_diff,ax);
    title(['插值节点数=',num2str(node_num)]);
    xlim([0,2*pi]);
    ylim([-1,1])
end
legend('sinx','插值',Location='bestoutside')

function cubic_Hermite_interpolation(x,y,y_diff,ax)
    node_num = length(x);
    for i = 1:node_num-1
        m = x(i);
        M = x(i+1);
        b = [y(i),y(i+1),y_diff(i),y_diff(i+1)]';
        A = [m^3,m^2,m^1,1;
            M^3,M^2,M^1,1;
            3*m^2,2*m,1,0;
            3*M^2,2*M,1,0
            ];
        coef = A\b;
        interpolation(coef,m,M,ax);
        hold on
    end
end

function interpolation(coef,m,M,ax)
    t = linspace(m,M,1000)';
    f = zeros(1000,1);
    for k = 1:1000 % 计算插值多项式在各点的函数值
        v = t(k,1);
        temp = ones(4,1);
        for p = 1:3
            temp(p) = v^(4-p);
        end
        f(k,1) = coef'*temp;
    end
    plot(ax,t,f,'b');
end