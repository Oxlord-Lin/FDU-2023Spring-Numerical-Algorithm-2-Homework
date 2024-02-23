clear;
clc;
runtime = zeros(2,1);
count = 0;
begin = 3;
last = 18;
for i = begin:last
    N = 2^i;
    count = count + 1;
    n1 = N/2; % u的最高次数
    n2 = n1;  % v的最高次数
    u = 10*randn(1,n1+1) + 10*1i*rand(1,n1+1); % 系数升幂排列
    v = 20*rand(1,n2+1) + 20*1i*randn(1,n2+1); % 系数升幂排列
    f = @() poly_product(u,v);
    t1 = timeit(f);
    runtime(1,count) = t1;
    g = @() conv(u,v);
    t2 = timeit(g);
    runtime(2,count) = t2;
end
%% 绘图
% my poly product function, log2(runtime/log2(n))~log2(n)
figure
h1 = plot(begin:last,log2(runtime(1,:)./(begin:last)),'b-o');
xlabel('log_2(n),n = the length of the vectors')
ylabel('log_2(runtime / log_2(n))')
legend('my poly product function')
saveas(gcf,'exercise-2-pic1.png')
% my poly product function, runtime/n~log2(n)
figure
h2 = plot(begin:last,runtime(1,:)./2.^(begin:last),'b--o');
xlabel('log_2(n), where n = the length of the vectors')
ylabel('runtime / n')
legend('my poly product function')
saveas(gcf,'exercise-2-pic2.png')
% conv log2(runtime/log2n)~log2(n)
figure
h3 = plot(begin:last,log2(runtime(2,:)./(begin:last)),'r-o');
xlabel('log_2(n),n = the length of the vectors')
ylabel('log_2(runtime / log_2(n))')
legend('conv')
saveas(gcf,'exercise-2-pic3.png')
%
% conv log2(runtime)~log2(n)
figure
h4 = plot(begin:last,log2(runtime(2,:)),'r-o');
xlabel('log_2(n),n = the length of the vectors')
ylabel('log2(runtime)')
legend('conv')
saveas(gcf,'exercise-2-pic4.png')


%% 
function y = poly_product(u,v)
% input: 按升幂排序的多项式系数向量，u和v
% output: 做乘法后的多项式系数向量y，升幂排序
    u = u(:);
    v = v(:);
    n1 = length(u)-1; % 最高次数n1
    n2 = length(v)-1; % 最高次数n2
    u = [u;zeros(n2,1)]; % 系数升幂排列
    v = [v;zeros(n1,1)]; 
    y = ifft(fft(u).*fft(v)); % 系数升幂排列
end