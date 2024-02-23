%%
for j = 1:4
    n = (j+1)*3; % 插值节点数
    fun = @(x) (1+25*(x.^2)).^(-1);
    x1 = linspace(-1,1,n)'; % 等间隔分布
    x2 = zeros(n,1);
    dominator = 2*n;
    for i = 1:n
        x2(i,1) = cos((2*i-1)*pi./dominator); % 切比雪夫插值节点
    end
    ax = nexttile;
    y1 = fun(x1);
    y2 = fun(x2);
    coef1 = ceof_Newton_basis(x1,y1);
    interpolation(ax,x1,coef1);
    hold on
    coef2 = ceof_Newton_basis(x2,y2);
    interpolation(ax,x2,coef2);
    t = linspace(-1,1,1000);
    f = fun(t);
    plot(ax,t,f,'--',LineWidth=1.5)
    % xlim([-1,1])
    legend('等距节点插值','Chebyshev节点插值','原函数图像')
    ylim([-1,1.5])
    title(['插值节点数为=',num2str(n)])

end
%%
function coef = ceof_Newton_basis(x,y)    
    y = y';
    x = x';
    n = length(x);
    divided_difference = zeros(n,n);
    divided_difference(1,:) = y; % 第一行作为0阶差商
    for i = 2:n
        for j = 1:n+1-i
            divided_difference(i,j) = (divided_difference(i-1,j+1) ...
                -divided_difference(i-1,j))/(x(j+i-1)-x(j));
        end
    end
    coef = divided_difference(:,1); % 牛顿插值的系数；
end

function interpolation(ax,x,coef)
% 下面是插值的部分
    n = length(x);
    m = min(x);
    M = max(x);
    t = linspace(m,M,1000)';
    f = zeros(1000,1);
    for k = 1:1000 % 计算插值多项式在各点的函数值
        v = t(k,1);
        temp = ones(n,1);
        for l = 2:n
            temp(l) = temp(l-1)*(v-x(l-1));
        end
        f(k,1) = coef'*temp;
    end
    plot(ax,t,f,LineWidth = 1.5);
end 