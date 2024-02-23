%% 通过二分法，找到临界值α
fun = @(x) -2*x + (1 + x^2) * atan(x); % 原函数
a = 1.25;
b = 1.5;
c = (a+b)/2;
while abs(fun(c)) > 1e-16 && abs(a-b)> 1e-16 % 二者之中只要有一个太小就会停止循环
    if fun(c)>0
        b = c;
        c = (a+b)/2;
    elseif fun(c) == 0
        break
    else 
        a = c;
        c = (a+b)/2;
    end
end
disp(['找到的阿尔法近似为：',num2str(c,15)])
disp(['迭代结束时|a-b|为：',num2str(abs(a-b))])