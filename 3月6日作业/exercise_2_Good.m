%% 用"Good" Broyden's method求解雨水的pH值
%% 准备部分
% 各平衡常数都抽象成常量
c1 = 1e-14;
c2 = 10^(-6.3 - 6 - 1.46) * 375;
c3 = 10^(-10.3);
% 用于记录误差变化情况
er1 = zeros(1,1); % 记录Good Broyden的收敛过程
er2 = zeros(1,1); % 记录Bad Broyden的收敛过程
% 初始状态
% 注意：初始状态要精心选取，否则会收敛到其他的解
% 根据化学常识可以判断，酸雨的pH大约在6附近，碳酸氢根的浓度约等于氢离子浓度
x = 1e-5;  % x是氢离子浓度
y = 1e-9;  % y是氢氧根浓度
z = 1e-5;  % z是碳酸氢根浓度
w = 1e-9;  % w是碳酸根浓度，根据化学常识，弱酸根离子的平衡常数一般相差1e3左右
v_0 = [x,y,z,w]';  % 将各离子浓度存储为列向量v
F = myfunc(v_0);
tiledlayout(2,1)
%% "Good" Broyden
v = v_0; 
inv_J = eye(4);
count = 0;
while norm(myfunc(v)) > 1e-16
    count = count + 1;
    delta_v = -inv_J*myfunc(v);
    delta_F = myfunc(v+delta_v) - myfunc(v);
    v = v + delta_v;
    inv_J = inv_J + ((delta_v - inv_J*delta_F)/(delta_v'*inv_J*delta_F))*delta_v'*inv_J;
    er1(count,1) = norm(myfunc(v));
end
pH = -log10(v(1,1));
disp(['雨水的pH=',num2str(pH)]);
ax = nexttile;
semilogy(ax,er1)
xlabel('迭代次数')
ylabel('误差（取对数）')
if pH>5 && pH<7 
    [~,~] = title('使用Good Broyden method误差随迭代次数的下降情况（正常情况）' , ['初始估计值为：','氢离子=',num2str(x),'，碳酸氢根=',num2str(z),'，碳酸根=',num2str(w),'，最终结果为pH=',num2str(-log10(v(1,1)))]);
else
    [~,~] = title('使用Good Broyden method误差随迭代次数的下降情况（异常情况）' , ['初始估计值为：','氢离子=',num2str(x),'，碳酸氢根=',num2str(z),'，碳酸根=',num2str(w),'，最终结果为pH=',num2str(-log10(v(1,1)))]);
end

%% "Bad" Broyden
% 效果并不好？
v = v_0; 
J = eye(4);
inv_J = eye(4);
count = 0;
while norm(myfunc(v)) > 1e-16
    count = count + 1;
    delta_v = - inv_J*myfunc(v);
    delta_F = myfunc(v+delta_v) - myfunc(v);
    v = v + delta_v;
    inv_J = inv_J + (delta_v - inv_J*delta_F)*delta_F'/norm(delta_F)^2;
    er2(count,1) = norm(myfunc(v));
end
pH = -log10(v(1,1));
disp(['雨水的pH=',num2str(pH)]);
ax = nexttile;
semilogy(ax,er2)
xlabel('迭代次数')
ylabel('误差（取对数）')
if pH>5 && pH<7 
    [~,~] = title('使用Bad Broyden method误差随迭代次数的下降情况（正常情况）' , ['初始估计值为：','氢离子=',num2str(x),'，碳酸氢根=',num2str(z),'，碳酸根=',num2str(w),'，最终结果为pH=',num2str(-log10(v(1,1)))]);
else
    [~,~] = title('使用Bad Broyden method误差随迭代次数的下降情况（异常情况）' , ['初始估计值为：','氢离子=',num2str(x),'，碳酸氢根=',num2str(z),'，碳酸根=',num2str(w),'，最终结果为pH=',num2str(-log10(v(1,1)))]);
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




