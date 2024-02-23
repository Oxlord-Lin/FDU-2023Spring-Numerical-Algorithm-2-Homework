% Radix-2 FFT
count = 1;
begin = 3;
last = 17;
run_time = zeros(last-begin+1,1); % 用于记录Radix-2 FFT的运行时间
for i = begin:last
    len = 2^i;
    vector = randn(len,1);
    % test
    % y = radix_2_fft(vector);
    % test ending
    f = @() radix_2_fft(vector);
    run_time(count) = timeit(f); % 计算调用radix_2_fft函数的时间
    count = count+1;
end
%% 绘图部分
h1 = plot(2.^(begin:last),run_time,'-o',LineWidth=1,Color='b'); 
hold on
% 绘制一条参考线
s = 2.^(begin:last);
s = s(:);
C = (s.*log2(s))\run_time;
h2 = plot(2.^(begin:last), C*s.*log2(s), '--', Color='r', LineWidth=1);
legend([h1,h2], 'runtime', ['y =',num2str(C),'*n*log_2n'], Location='southeast');
xlabel('the length of vector');
title('the runtime of Radix-2 FFT against the length of vector')
%% 另一种绘图方式，更便于看出nlogn
xx = begin:last; % 取对数,log2(n)
yy = run_time'./xx; % runtime./log2(n) ~ O(n)
h3 = plot(xx,yy,'-o');
xlabel('log2(n), where n = the length of vectors')
ylabel('runtime / log2(n) ~ O(n)');

%%
function x = radix_2_fft(x)
    x = x(:);
    n = length(x);
    % 先做排序，排序的复杂度是O(n)
    for i = log2(n): -1 : 2 % i是每"段"的长度取对数
        block = 2^i; % block是每一"段"的长度
        for j = 1 : n/block % 第j段
            a = (j-1)*block;
            b = j*block;
            u = x(a + 1 : 2: b-1); % 奇数部分
            v = x(a + 2 : 2: b); % 偶数部分
            x(a + 1 : b) = [u;v];% 奇数部分放在前，偶数部分放在后
        end
    end
    % 再做剩下的乘法部分
    for i = 1 : log2(n) % i是每"段"的长度取对数
        block = 2^i; % block是每一"段"的长度
        for j = 1 : n/block % 第j段
            omega = exp(-2*pi*1i/block);
            k = [0:block/2-1]';
            w = omega.^k;
            a = (j-1)*block;
            b = j*block;
            u = x(a + 1 : a+block/2); % 奇数部分
            v = w .* x(a + block/2 + 1 : b); % 偶数部分，做Hadamard积
            x(a + 1 : b) = [u+v;u-v]; % 这样就完成了一次局部的FFT
        end
    end
end













