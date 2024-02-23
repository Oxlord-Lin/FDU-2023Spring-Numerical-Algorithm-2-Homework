%% Gaussian smoother
clear;
clc;
a = 2;
N = 10000; % 切成N份
xx = linspace(-a,a,N+1);
xx = xx(1:N);
f = myfunc1(xx);
miu = 0;
% variance = 0.05
variance1 = 0.05;
g = gauss(xx,miu,variance1);
h1 = conv(f,g);
h1 = h1(N/2+2:N+N/2+1); % 【注意，这里的截断颇有讲究】广义离散卷积，取N项
h1 = 2*a/N*h1; % 积分离散化时有一个系数
% variance = 0.01
variance2 = 0.01;
g = gauss(xx,miu,variance2);
h2 = conv(f,g);
h2 = h2(N/2+2:N+N/2+1); % 【注意，这里的截断颇有讲究】广义离散卷积，取N项
h2 = 2*a/N*h2; % 积分离散化时有一个系数%
% variance = 0.001
variance3 = 0.001;
g = gauss(xx,miu,variance3);
h3 = conv(f,g);
h3 = h3(N/2+2:N+N/2+1); % 【注意，这里的截断颇有讲究】广义离散卷积，取N项
h3 = 2*a/N*h3; % 积分离散化时有一个系数
%% 绘图部分
figure
p = bar(xx,f,FaceColor=[0.95,0.95,0]);
hold on
p1 = plot(xx,h1,'m',LineWidth=1.5,LineStyle=':');
p2 = plot(xx,h2,'r',LineWidth=1.5,LineStyle='--');
p3 = plot(xx,h3,'c',LineWidth=1.5,LineStyle='-');
legend([p,p1,p2,p3],'原本的非连续函数',['\sigma^2=',num2str(variance1)],['\sigma^2=',num2str(variance2)],['\sigma^2=',num2str(variance3)],Location='best')
title(['高斯磨光算子效果展示'])
saveas(gcf,'exercise-3-高斯磨光.png')




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

function g = gauss(x,mean,variance)
    n = length(x);
    g = zeros(n,1);
    for i = 1:n
        g(i) = (1/sqrt(2*pi*variance))*exp(-(x(i)-mean)^2/(2*variance));
    end
end
