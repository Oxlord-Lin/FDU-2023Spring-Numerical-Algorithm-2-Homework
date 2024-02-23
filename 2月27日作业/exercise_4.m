%% 用牛顿法求解初始估计值x_0 = 3时的方程 1+cosx=0
% 已知最靠近x0 = 3的零点是pi
fun = @(x) 1 + cos(x);
x = 3;
er = zeros(1,30);
count = 1;
er(1,count) = abs(x-pi);
while fun(x) > 1e-15
    count = count + 1;
    x = x + fun(x)/sin(x);
    er(1,count) = abs(x-pi);
end
% 对误差变化作半对数图
semilogy(er(1,1:count));
title('误差（取对数）随迭代次数变化的情况')
xlabel('迭代次数')
ylabel('log_{10}|x_k-x_*|')

% 观察线性收敛时的比例系数
relative_er = zeros(1,count-1);
for j = 1:count-1
    relative_er(1,j) = er(1,j+1)/er(1,j);
end
figure
plot(relative_er);
xlabel('迭代次数')
ylabel('$$\frac{|x_{k+1}-x_*|} {|x_k-x_*|}$$',Interpreter='latex') 
% 比例系数大约为0.5，具体的证明可参看实验报告