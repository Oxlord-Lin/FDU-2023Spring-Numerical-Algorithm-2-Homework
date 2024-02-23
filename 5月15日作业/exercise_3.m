h = 0.5; % 初始步长
n = 55; % 外插的阶数
x = 1; % 求何处的导数
inter_results = Richardson_extrapolation(x,h,n);
%% 绘图
plot(1:n,inter_results,'b->',LineWidth=1)
title('f(x)=x^3e^x在x=1处的导数计算值随Richardson外插阶数的变化')
xlabel('外插阶数')
ylabel('导数计算值')
legend('导数计算值',Location='best')
saveas(gcf,'exercise-3 导数计算值.png')
figure
real_dev = exp(x)*(3*x^2 + x^3);
semilogy(1:n,abs(inter_results - real_dev),'r--*',LineWidth=1);
title('f(x)=x^3e^x在x=1处的导数计算值的误差随Richardson外插阶数的变化')
xlabel('外插阶数')
ylabel('导数计算的误差值（取对数）')
legend('导数计算值的误差',Location='best')
saveas(gcf,'exercise-3 导数计算误差值.png')

function inter_results = Richardson_extrapolation(x,h,n)
    if n == 1
        inter_results = (func(x+h) - func(x-h))/ (2*h);
    elseif n > 1
        Table = zeros(n,n);
        for i = 1:n
            Table(i,1) = (func( x + h/2^(i-1) ) - func( x - h/2^(i-1)) ) / (2 * h/2^(i-1) );
        end
        for j = 2:n
            for i = 1:n-j+1
                Table(i,j) = ( 2^(j-1)*Table(i+1,j-1) - Table(i,j-1) ) / (2^(j-1) - 1);
            end
        end
        inter_results = Table(1,:);
        return
    else
        error('n小于1，错误！');
    end
end

function f = func(x)
    f = x^3 * exp(x);
end