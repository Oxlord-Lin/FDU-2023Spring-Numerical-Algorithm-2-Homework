% Friedrichs smoother
clear;
clc;
a = 2;
N = 10000; % 切成N份
xx = linspace(-a,a,N+1);
xx = xx(1:N);
f = myfunc1(xx);
miu = 0;
% epsilon = 0.4
epsilon1 = 0.4;
g1 = Friedrichs_molliferis(xx,epsilon1);
h1 = conv(f,g1);
h1 = h1(N/2+2:N+N/2+1); % 【注意，这里的截断颇有讲究】广义离散卷积，取N项
h1 = 2*a/N*h1; % 积分离散化时有一个系数
% epsilon = 0.2
epsilon2 = 0.2;
g2 = Friedrichs_molliferis(xx,epsilon2);
h2 = conv(f,g2);
h2 = h2(N/2+2:N+N/2+1); % 【注意，这里的截断颇有讲究】广义离散卷积，取N项
h2 = 2*a/N*h2; % 积分离散化时有一个系数
% epsilon = 0.1
epsilon3 = 0.1;
g3 = Friedrichs_molliferis(xx,epsilon3);
h3 = conv(f,g3);
h3 = h3(N/2+2:N+N/2+1); % 【注意，这里的截断颇有讲究】广义离散卷积，取N项
h3 = 2*a/N*h3; % 积分离散化时有一个系数


%% 绘图部分
figure
p = bar(xx,f,FaceColor=[0.95,0.95,0]);
hold on
p1 = plot(xx,h1,'m',LineWidth=1.5,LineStyle=':');
p2 = plot(xx,h2,'r',LineWidth=1.5,LineStyle='--');
p3 = plot(xx,h3,'c',LineWidth=1.5,LineStyle='-');
legend([p,p1,p2,p3],'原本的非连续函数',['\epsilon=',num2str(epsilon1)],['\epsilon=',num2str(epsilon2)],['\epsilon=',num2str(epsilon3)],Location='best')
title(['Friedrichs磨光算子效果展示'])
saveas(gcf,'exercise-3-Friedrichs磨光.png')

%% 函数
function f = myfunc1(x)
    n = length(x);
    f = zeros(n,1);
    for i = 1:n
        temp = x(i);
        if -2<temp && temp<-1
            f(i) = -0.5*exp(sin(-3*x(i)))*sign(cos(14*cos(x(i))));
%               f(i) = -1;
        elseif -1<temp && temp<0
            f(i) = 0.3*sign(cos(10*x(i)));
%               f(i) = -1;
        elseif 0<temp && temp<1
            f(i) = -0.3*x(i).^2*sign(3*cos(10*x(i)));
%             f(i) = 1;
        elseif 1<temp && temp<2
            f(i) = sin(2*x(i))*log(x(i))*sign(cos(exp(x(i))*cos(x(i))));
%             f(i) = 1;
        else 
            f(i) = 0;
        end
    end
    f = real(f);
end

function g = Friedrichs_molliferis(x,epsilon)
    C = 0.4439938161680794378;
    n = length(x);
    g = zeros(n,1);
    for i = 1:n
        if abs(x(i))<= epsilon
            g(i) = (1/C)*(1/epsilon)*exp(1/( x(i)^2 / epsilon^2 - 1)); % 这里很危险，容易出inf
        else
            g(i) = 0;
        end
        if isinf(g(i)) % 一定要有这步，防止inf
            g(i) = 0;
        end
    end
end
