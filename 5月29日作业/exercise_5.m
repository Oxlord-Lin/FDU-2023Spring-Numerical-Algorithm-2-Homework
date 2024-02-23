clear;
clc;
e = exp(1);
alpha = -2 - (-2*e + 2) / (e - 1/e);
beta = (-2*e + 2) / (e - 1/e);
u = @(x) x.^2 + 2 + alpha*exp(x) + beta*exp(-x);
max_error_FDM = [];
for j = 1:9
    n = 2^j;
    xx = linspace(0,1,n+1);
    u_exact = u(xx);
    u_FDM = FDM(0,1,n);
    max_error_FDM = [max_error_FDM, max(abs(u_exact(:)-u_FDM(:)))];
    ax = nexttile;
    h1 = plot(ax,xx,u_exact,'r-',LineWidth=1);hold on
    h2 = plot(ax,xx,u_FDM,'b--',LineWidth=1);
    title(['FDM, 节点个数=',num2str(n+1)]);
end
legend([h1,h2],'精确值','估计值',Location='bestoutside');
figure
loglog(2.^(1:9),max_error_FDM,'r-o',LineWidth=1);
title('FDM误差随节点数n的变化（双对数图），误差为O(n^{-2})');
xlabel('log(n)')
xlabel('log(error)')
saveas(gcf,'exercise-5-FDM-error.png')
 
%% 
max_error_FEM = [];
for j = 1:9
    n = 2^j;
    xx = linspace(0,1,n+1);
    u_exact = u(xx);
    u_FEM = FEM(0,1,n-1);
    max_error_FEM = [max_error_FEM, max(abs(u_exact(:)-u_FEM(:)))];
    ax = nexttile;
    h1 = plot(ax,xx,u_exact,'r-',LineWidth=1);hold on
    h2 = plot(ax,xx,u_FEM,'b--',LineWidth=1);
    title(['FEM, 节点个数=',num2str(n+1)]);
end
legend([h1,h2],'精确值','估计值',Location='bestoutside');
figure
loglog(2.^(1:9),max_error_FEM,'r-o',LineWidth=1);
title('FEM误差随节点数n的变化（双对数图），误差为O(n^{-2})');
xlabel('log(n)')
xlabel('log(error)')
saveas(gcf,'exercise-5-FEM-error.png')



%% 
function U = FDM(u_start,u_final,n)
% input:u_start:左端初值，u_fianl：右端初值，n：最终会把区间切割成n段，一共有n+1个均分点
% output: U：n+1维向量，在n+1个均分点上的函数估计值
    x = linspace(0,1,n+1);
    x = x(:);
    h = 1/n;
    e = ones(n-1,1);
    A = spdiags([-e (2+h^2)*e -e],-1:1,n-1,n-1);
    b = h^2 * (x(2:n)).^2;
    b(1) = b(1) + u_start;
    b(end) = b(end) + u_final;
    u = A\b;
    U = [u_start;u;u_final];
end


function C = FEM(u_start,u_final,n)
% input:u_start:左端初值，u_fianl：右端初值，n:最终会把区间切割成n+1段，一共有n+2个均分点
% output: C：n+2维向量，在n+2个均分点上的函数估计值
    C = zeros(n+2,1);
    C(1) = u_start;
    C(end) = u_final;
%     x = linspace(0,1,n+2);
    h = 1/(n+1);
    e = ones(n,1);
    A = spdiags([(h/6 - 1/h)*e, (2/3*h + 2/h)*e, (h/6 - 1/h)*e],-1:1,n,n);
    b = fun2(n); % 自己写了一个函数，因为有限元和x.^2的积分太复杂了
    b(1) = b(1) - (h/6 - 1/h)*u_start;
    b(end) = b(end) - (h/6 - 1/h)*u_final;
    c = A\b;
    C(2:end-1) = c;
end

function b = fun1(n)
    x = linspace(0,1,n+2);
    x = x(:);
    x1 = x(1:end-2);
    x2 = x(2:end-1);
    x3 = x(3:end);
    % 记得写系数
    b = (n+1)^(-3)*(1/4*(2*x2.^4 - x1.^4 - x3.^4) + 1/3*(0:n-1)'.*(-x2.^3 + x1.^3) + 1/3*(2:n+1)'.*(x3.^3 - x2.^3));
end

function b = fun2(n)
    b = zeros(n,1);
    for i = 1:n
        b(i) = (1/4)*(n+1)^(-3)*(2*i^4 - (i-1)^4 - (i+1)^4) + (i-1)*(1/3)*(n+1)^(-3)*(-i^3+(i-1)^3) +  (i+1)*(1/3)*(n+1)^(-3)*((i+1)^3 - (i)^3);
    end 
end