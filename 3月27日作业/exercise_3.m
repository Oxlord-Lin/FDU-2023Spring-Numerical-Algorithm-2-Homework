fun = @(x) x^64 - 0.1; % 原函数

%% regula falsi方法求根
a = 0;
y_a = fun(a);
b = 1;
y_b = fun(b);
c = (b*y_a - a*y_b)/(y_a - y_b);
regula_falsi_abs_er = zeros(1,1); % 用于记录收敛过程
regula_falsi_abs_er(1,1) = abs(fun(c));
count = 1;
while abs(fun(c))> 1e-12 % 使用abs(fun(c))作为误差判断依据
    count = count + 1;
    f_c = fun(c);
    if f_c == 0 
        break
    elseif f_c > 0
        b = c;
        f_b = f_c;
    else
        a = c;
        f_a = f_c;
    end
    c = (b*y_a - a*y_b)/(y_a - y_b);
    regula_falsi_abs_er(count,1) = abs(fun(c)); % 记录误差
end
semilogy(regula_falsi_abs_er,':',LineWidth=1,Color='g')
xlabel("迭代次数")
ylabel('误差（取对数）')
hold on
% 选取残差曲线上的下述诸点，用最小二乘法进行直线拟合
nodes = [[106,6.08887940067859e-13]
[105,2.78471690151605e-12]
[102,1.89699367325602e-11]
[101,1.08730247028177e-11]
[98,1.21812546161060e-10]
[91,5.68332797490356e-10]
[90,2.87557963551954e-10]
[87,2.89414797682852e-09]
[82,9.04808289403647e-09]
[77,1.33148250605619e-07]
[67,6.99226859510183e-07]
[66,1.69897529780916e-07]
[62,4.38560398458021e-06]
[57,8.86155040026593e-06]
[54,6.61731602582538e-05]
[46,0.000431707413325017]
[45,0.000210983053052935]
[42,0.00215774526799048]
[33,0.0149187083131281]
[32,0.00686672734992283]
[29,0.0543972129900900]
[16,0.0999979872725661]
[1,0.100000000000000]];

x = nodes(:,1);
y = nodes(:,2);
plot(x,y,'*',Color='r')
y = log10(y); % 直接用原始数据会出现条件数很大的情况，因此先取对数，再还原回去
n = length(y);
A = [x,ones(n,1)];
t = A\y;
k = t(1); % 斜率
b = t(2); % 截距
xx = 0:0.01:max(x);
yy = k*xx + b;
for i = 1:length(yy)
    yy(i) = 10^(yy(i));
end
plot(xx,yy,LineWidth=1,Color='b')
legend('the history of residual using regula falsi','选取的点','最小二乘法拟合的指数类曲线（图中取了对数）')
[~,~] = title(['the history of residual using regula falsi and the curve that fits the history'], ...
    ['k=',num2str(k),' b=',num2str(b)]);




%% 二分法求根
a = 0;
b = 1;
c = 0.5*(a+b);
bisection_abs_er = zeros(1,1);
bisection_abs_er(1,1) = abs(fun(c));
count2 = 1;
while abs(fun(c)) > 1e-12
    count2 = count2 + 1;
    if fun(c)>0
        b = c;
        c = (a+b)/2;
    elseif fun(c) == 0
        break
    else 
        a = c;
        c = (a+b)/2;
    end
    bisection_abs_er(count2,1) = abs(fun(c));
end
figure
semilogy(bisection_abs_er,'--',LineWidth=1,Color='g')
xlabel("迭代次数")
ylabel('误差（取对数）')
hold on
% 选取残差曲线上的下述诸点，用最小二乘法进行直线拟合
nodes = [[17,4.41553253232863e-05]
[41,5.77218828290427e-13]
[40,3.59422214213367e-12]
[37,9.62821489203236e-12]
[38,1.45077838631380e-11]
[35,1.54444165789691e-10]
[34,3.86437687405206e-11]
[32,3.47532086442115e-10]
[30,1.89223550717266e-09]
[18,1.88390574820674e-05]
[28,2.04286747101934e-08]
[27,4.28657831452739e-09]
[23,6.46885267421848e-07]
[25,5.37171023912153e-08]
[21,3.01959256282514e-06]
[20,1.44004849628177e-07]
[19,6.18328849505412e-06]
[15,0.000196185449859360]
[14,6.47090412728146e-06]
[13,0.000399246096061054]
[12,0.000410573902691724]
[9,0.0110951750290665]
[8,0.00121554750045422]
[7,0.0219291200033380]
[6,0.0536989216621620]
[5,0.0310840324784751]
[1,0.100000000000000]];

x = nodes(:,1);
y = nodes(:,2);
plot(x,y,'*',Color='r')
y = log10(y); % 直接用原始数据会出现条件数很大的情况，因此先取对数，再还原回去
n = length(y);
A = [x,ones(n,1)];
t = A\y;
k = t(1); % 斜率
b = t(2); % 截距
xx = 0:0.01:max(x);
yy = k*xx + b;
for i = 1:length(yy)
    yy(i) = 10^(yy(i));
end
plot(xx,yy,LineWidth=1,Color='b')
legend('the history of residual using bisection','选取的点','最小二乘法拟合的指数类曲线（图中取了对数）')
[~,~] = title(['the history of residual using bisection and the curve that fits the history'], ...
    ['k=',num2str(k),' b=',num2str(b)]);