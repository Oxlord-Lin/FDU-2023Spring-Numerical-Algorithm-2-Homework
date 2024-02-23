x = linspace(-10,10,4000)';
y1 = zeros(4000,1);
y2 = y1;
N = 10;
for j = 1:4000
    y1(j,1) = my_sin(x(j,1),N);
end
for j = 1:4000
    t = x(j,1);
    k = round(t/pi);
    t = t - k*pi;
    y2(j,1) = my_sin(t,N);
    if mod(k,2) == 1
        y2(j,1) = y2(j,1)*(-1);
    end

end
figure
plot(x,y1);
hold on;
plot(x,y2);
legend('方法a','方法b')
[~,~] = title('用方法a和方法b得到的sinx图像',['阶数：',num2str(N)]);

figure
semilogy(x,abs(y1-y2));




function result = my_sin(theta,N)
result = 0;
for n =0:N
    result = result +  (-1)^n * (theta)^(2*n+1)/factorial(2*n+1);
end
end