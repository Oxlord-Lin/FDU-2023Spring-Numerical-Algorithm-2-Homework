%% 用多元牛顿法求解雨水的pH值
% 各平衡常数都抽象成常量
c1 = 1e-14;
c2 = 10^(-6.3 - 6 - 1.46) * 375;
c3 = 10^(-10.3);
% 用于记录误差变化情况
er = zeros(1,1);
% 初始状态
x = 1e-6;  % x是氢离子浓度
y = 1e-6;  % y是氢氧根浓度
z = 1e-6;  % z是碳酸氢根浓度
w = 1e-6;  % w是碳酸根浓度
v = [x,y,z,w]';  % 将各离子浓度存储为列向量v
F = myfunc(v);
v_star = 1.0e-05 *[0.255480112575731;
                0.000391419899544;
                0.255078684681510;
                0.000005003997338]; % v_star 已经提前存储好（在上一轮求得），可以认为是精确解
er(1,1) = norm(v-v_star); 
% 牛顿法
count = 1;
while norm(F) > 1e-16
    count = count + 1;
    F = myfunc(v);
    F_d = myfunc2(v);
    step = F_d\myfunc(v);
    v = v - step;
    er(count,1) = norm(v-v_star);
end
pH = -log10(v(1,1));
disp(['雨水的pH值为：',num2str(pH)]);
semilogy(er);
xlabel('迭代次数')
ylabel('误差（取对数）')
if pH>5 && pH<7 
    [~,~] = title('使用Newton method误差随迭代次数的下降情况（正常情况）' , ['初始估计值为：','氢离子=',num2str(x),'，氢氧根离子=',num2str(y),'，碳酸氢根=',num2str(z),'，碳酸根=',num2str(w),'，最终结果为pH=',num2str(-log10(v(1,1)))]);
else
    [~,~] = title('使用Newton method误差随迭代次数的下降情况（异常情况）' , ['初始估计值为：','氢离子=',num2str(x),'，氢氧根离子=',num2str(y),'，碳酸氢根=',num2str(z),'，碳酸根=',num2str(w),'，最终结果为pH=',num2str(-log10(v(1,1)))]);
end
%% 计算F
function F = myfunc(v)
c1 = 1e-14;
c2 = 10^(-6.3 - 6 - 1.46) * 375;
c3 = 10^(-10.3);
x = v(1,1);
y = v(2,1);
z = v(3,1);
w = v(4,1);
F = [x*y - c1;
    x*z - c2;
    x*w/z - c3;
    x-y-z-2*w;];
end
%% 计算F的导数
function F_d = myfunc2(v)
x = v(1,1);
y = v(2,1);
z = v(3,1);
w = v(4,1);
F_d =  [y, x, 0, 0;
        z, 0, x, 0;
        w/z, 0, -x*w/z^2, x/z;
        1, -1, -1, -2;];
end