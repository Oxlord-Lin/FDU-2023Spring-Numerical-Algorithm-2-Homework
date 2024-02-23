%% 绘制不同插值节点的图像与sinx的对比图 
for n = 3:6
    x = linspace(0,2*pi,n);
    y = sin(x);
    ax = nexttile;
    interpolation(x,y); % v中存储的是多项式的系数
    hold on
    t = linspace(0,2*pi,1000);
    f = sin(t);
    plot(t,f,'--',LineWidth=1.5)
    xlim([0,2*pi])
    ylim([-1.2,1.2])
    legend(['n=',num2str(n)],'sin(x)')
    title(ax,['对sin(x)进行插值，节点数为=',num2str(n)])
end
%% 绘制差异图
figure
for n = 3:6
    x = linspace(0,2*pi,n);
    y = sin(x);
    ax = nexttile;
    interpolation2(x,y); % v中存储的是多项式的系数
    xlim([0,2*pi])
    ylim([0,1])
    title(ax,['对sin(x)进行插值的误差，节点数=',num2str(n)])
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
    % 下面是绘图，统一取[0,2*pi]之间的1000个点
    t = linspace(0,2*pi,1000);
    f = zeros(1000,1);
    for j = 1:1000
        f(j,1) = myfunc(v,t(1,j)); % f是处于t处的多项式函数值
    end
    plot(t,f,LineWidth=1.5);
    hold on
end

function interpolation2(x,y) % 用于绘制差异图
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
    % 下面是绘图，统一取[0,2*pi]之间的1000个点
    t = linspace(0,2*pi,1000);
    f = zeros(1000,1);
    for j = 1:1000
        f(j,1) = abs(myfunc(v,t(1,j))-sin(t(1,j))); % f是处于t(1,j)处的多项式函数值
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
