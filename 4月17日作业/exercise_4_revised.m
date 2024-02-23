n = 2; % 一共n+1个控制点，生成n阶Bezier曲线
points = randn(n+1,2); % 用户可以根据需要生成任意多个随机控制点，每一行表示一个点的横、纵坐标。

Bezier(points);

function Bezier(points)
    x = points(:,1);
    y = points(:,2);
    n = length(x) - 1; % 共n+1个点，生成n阶Bezier曲线
    C = zeros(n+1,1); % C用于存储二项式系数
    for i = 0:n
        C(i+1) = nchoosek(n,i);
    end
    temp_x = x.*C;
    temp_y = y.*C;
    tt = linspace(0,1,2000);
    xx = tt;
    yy = tt;
    for i = 1:2000
        temp_t = vector_t(tt(i),n); % 得到一个关于tt(i)的n+1维向量
        xx(i) = temp_t' * temp_x;
        yy(i) = temp_t' * temp_y;
    end
    plot(xx,yy,LineStyle="-",LineWidth=1); % 绘制Bezier曲线
    hold on
    scatter(x,y,'filled'); % 绘制控制点
    legend('Bezier曲线','控制点')
    title([num2str(n),'阶Bezier曲线，',num2str(n+1),'个控制点'])
end

function vector = vector_t(t,n)
    vector = zeros(n+1,1);
    for i = 1:n+1
        vector(i) = t^(i-1) * (1-t)^(n-i+1);
    end
end