%% 绘制对比图
for n = 1:6
    space = 3*n;
    x = linspace(-1,1,space);
    fun = @(x) 1./(1+25*x.^2);
    y = fun(x);
    ax = nexttile;
    interpolation(x,y); % v中存储的是多项式的系数
    hold on
    t = linspace(-1,1,1000);
    f = fun(t);
    plot(t,f,'--',LineWidth=1.5)
    ylim([-2,2])
%     legend(['n=',num2str(space)],'f(x)=(1+25x^2)^{-1}')
    title(ax,['节点数为=',num2str(space)])
end


%% 绘制差异图
figure
for n = 1:6
    space = 3*n;
    x = linspace(-1,1,space);
    fun = @(x) 1./(1+25*x.^2);
    y = fun(x);
    ax = nexttile;
    interpolation2(x,y); % v中存储的是多项式的系数
    ylim([0,5])
    title(ax,['对函数f(x)=(1+25x^2)^{-1}插值误差，','节点数=',num2str(space)])
end

%%
function interpolation(x,y)
    n = length(x);
    v = zeros(n,1); % 用于存储多项式系数
    for i = 1:n
        x_part = [x(1:i-1),x(i+1:end)]';
        dom =  1;% dom是分母
        for j = 1:n-1
            dom = dom*(x(i)-x_part(j));
        end
        v = v + poly(x_part)'.*y(i)./dom; %poly返回首一多项式系数
    end
    % 下面是绘图，统一取[-1,1]之间的1000个点
    t = linspace(-1,1,1000);
    f = zeros(1000,1);
    for j = 1:1000
        f(j,1) = myfunc(v,t(1,j)); % f是处于t处的多项式函数值
    end
    plot(t,f,LineWidth=1.5);
    hold on
end

function interpolation2(x,y)
    fun = @(x) 1./(1+25*x.^2);
    n = length(x);
    v = zeros(n,1); % 用于存储多项式系数
    for i = 1:n
        x_part = [x(1:i-1),x(i+1:end)]';
        dom =  1;% dom是分母
        for j = 1:n-1
            dom = dom*(x(i)-x_part(j));
        end
        v = v + poly(x_part)'.*y(i)./dom; %poly返回首一多项式系数
    end
    % 下面是绘图，统一取[-1,1]之间的1000个点
    t = linspace(-1,1,1000);
    f = zeros(1000,1);
    for j = 1:1000
        f(j,1) = abs(myfunc(v,t(1,j))-fun(t(1,j))); % f是处于t处的多项式函数值
    end
    plot(t,f,LineWidth=1.5);
    hold on
end


function f = myfunc(v,t)
    n = length(v);
    temp = ones(n,1);
    for i = 1:n
        temp(i) = t^(n-i);
    end
    f = v'*temp;
end
