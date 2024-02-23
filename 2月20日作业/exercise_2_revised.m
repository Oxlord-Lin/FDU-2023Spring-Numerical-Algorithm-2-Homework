%%
num = 1000;
x = linspace(-10,10,num)';
y1 = zeros(num,1);
y2 = y1;
N1 = zeros(num,1); %用来观察截断项数的变化情况
N2 = zeros(num,1); 
for j = 1:num
    [y1(j,1),N1(j,1)] = my_sin(x(j,1));
end
for j = 1:num
    t = x(j,1);
    k = round(t/pi);
    t = t - k*pi;
    [y2(j,1),N2(j,1)] = my_sin(t);
    if mod(k,2) == 1
        y2(j,1) = y2(j,1)*(-1);
    end

end

%%
figure
plot(x,y1,'--',LineWidth=1.5);
hold on;
plot(x,y2,'-');
legend('方法a','方法b')
title('用方法a和方法b得到的sinx图像');

figure
semilogy(x,abs(y1-y2));
title('两种方法的差距（取对数）');
xlabel('x');

figure
plot(x,N1);
title('N1随着x的变化情况')

figure
plot(x,N2);
title('N2随着x的变化情况')



function [result,N] = my_sin(theta)
% 先确定N
N = ceil(abs(theta)/2);
u = 1e-16;
RN = (abs(theta))^(2*N+3)/factorial(2*N+3);
gama_3N = 3*N*u/(1-3*N*u);
for t = 1:100
    if RN < gama_3N*(1 + RN)
        break
    else
        N = N+10;
        RN = (abs(theta))^(2*N+3)/factorial(2*N+3);
        gama_3N = 3*N*u/(1-3*N*u);
    end
end
result = 0;
for n =0:N
    result = result +  (-1)^n * (theta)^(2*n+1)/factorial(2*n+1);
end
end