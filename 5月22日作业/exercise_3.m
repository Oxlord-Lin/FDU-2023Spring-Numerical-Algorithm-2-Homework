%% 计算 f(x) = exp(x) 在[0，1]上的积分
n = 10;
a = 0;
b = 1;
r1 = romberg(@exp, 0, 1, n);
h = (b-a)./(2.^(0:n-1)); % 所有用到过的步长
exact_result_1 = exp(1) - 1;
er_1 = abs(r1 - exact_result_1);

% 第1列是复合梯形公式，误差是 O（h^(2*1)）线性收敛
h1 = semilogy(er_1(:,1),':o',LineWidth=1); hold on

% 第2列是复合Simpson公式，误差是 O（h^(2*2)）线性收敛
h2 = semilogy(er_1(2:end,2),'--s',LineWidth=1);

% 第3列是复合Cotes公式，误差是 O（h^(2*3)）线性收敛
h3 = semilogy(er_1(3:end,3),'->',LineWidth=1);

% Romberg积分的最优结果是每行最右边的元素，误差是 O（ ( (b-a)/2^(j-1) )^(2*j) ）
% 在本题中也即O（ (1/2)^(2*j*(j-1)) ），也即二次收敛
h4 = semilogy(diag(er_1),'-.h',LineWidth=1);
xlabel('iteration',Interpreter='latex',FontSize=14)
ylabel('$log(error)$',Interpreter='latex',FontSize=14)
title('四种积分对f(x)=exp(x)在[0,1]上积分的收敛情况')
legend([h1,h2,h3,h4],'复化梯形公式','复化Simpson公式','复化Cotes公式','Romberg积分')
saveas(gcf,'exercise3-exp(x).png')

%% 计算 f(x) = x^(1.5) 在[0，1]上的积分
n = 20;
a = 0;
b = 1;
f = @(x) x^1.5;
r2 = romberg(f,a,b,n);
h = (b-a)./(2.^(0:n-1)); % 所有用到过的步长
exact_result_2 = 1/2.5;
er_2 = abs(r2 - exact_result_2);

% 第1列是复合梯形公式，误差是 O（h^(2*1)）线性收敛
figure
h1 = semilogy(er_2(:,1),':o',LineWidth=1); hold on

% 第2列是复合Simpson公式，误差是 O（h^(2*2)）线性收敛
h2 = semilogy(er_2(2:end,2),'--s',LineWidth=1);

% 第3列是复合Cotes公式，误差是 O（h^(2*3)）线性收敛
h3 = semilogy(er_2(3:end,3),'->',LineWidth=1);

% Romberg积分的最优结果是每行最右边的元素，误差是 O（ ( (b-a)/2^(j-1) )^(2*j) ）
% 在本题中也即O（ (1/2)^(2*j*(j-1)) ），也即二次收敛
h4 = semilogy(diag(er_2),'-.h',LineWidth=1);
xlabel('iteration',Interpreter='latex',FontSize=14)
ylabel('$log(error)$',Interpreter='latex',FontSize=14)
title('四种积分对f(x)=x^{1.5}在[0,1]上积分的收敛情况')
legend([h1,h2,h3,h4],'复化梯形公式','复化Simpson公式','复化Cotes公式','Romberg积分')
saveas(gcf,'exercise-3-x^{1.5}.png')





function r = romberg(f,a,b,n)
% input： f：函数句柄, [a,b]：积分区间 n：Romberg表的行数
% ouput： n行Romberg表
    h = (b-a)./(2.^(0:n-1)); % 所有要用到的步长都在这里
    r(1,1) = (b-a) * (f(a) + f(b))/2;
    for i = 2:n
        subtotal = 0;
        for k = 1:2^(i-2)
            subtotal = subtotal + f(a + (2*k - 1)*h(i));
        end
        r(i,1) = r(i-1,1)/2 + h(i)*subtotal; % 算出第一列
        for j = 2:i % 往右递推
            r(i,j) = (4^(j-1)*r(i,j-1) - r(i-1,j-1)) / (4^(j-1) - 1);
        end
    end
end