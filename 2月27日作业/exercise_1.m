%% verify
y1 = linspace(-1,10,1000);
myfunc = @(y) y.*exp(y);
x1 = myfunc(y1);

%% 使用牛顿法找到每个x对应的y（逐点反解）
x2 = linspace(-1*exp(-1),10*exp(10),1000);
y2 = zeros(1,1000);
for i = 1:1000
    y2(1,i) = root_finder_1(x2(1,i));
end

%% 画图
plot(x1,y1,'--',LineWidth=1.5); % verify
hold on
plot(x2, y2,':',LineWidth=1); % 逐点牛顿法反解
legend('准确的函数图像','通过牛顿法逐点反解得到的函数图像');

%% 用牛顿法求根
function root = root_finder_1(c)
    count = 0;
    fun = @(y) y*exp(y) - c;
    fun_de = @(y) (y+1)*exp(y);
    root = 10;
    res = fun(root);
    while abs(res) > 1e-15  && count < 200 % 这里加了一个迭代次数保护
        root = root - fun(root)/fun_de(root);
        res = fun(root);
        count = count + 1;
    end
end
