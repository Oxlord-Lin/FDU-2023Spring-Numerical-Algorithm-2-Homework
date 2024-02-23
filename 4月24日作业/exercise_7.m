% Radix-3 FFT
count = 1;
begin = 2;
last = 13;
runtime = zeros(last-begin+1,1);
for k = begin:last
    n = 3^k;
    vector = randn(n,1);
    % test
    % y = radix_3_fft(vector);
    % test ending
    tic;
    radix_3_fft(vector);
    runtime(count) = toc;
    count = count+1;
end
%% 绘图
xx = 3.^(begin:last);
h1 = plot(xx,runtime,'-o',LineWidth=1); 
hold on
xlabel('the length of the vector')
ylabel('runtime')
title('the runtime of my Radix-3 FFT against the length of the vector')
% 绘制辅助参考线
% 绘制一条参考线
s = 3.^(begin:last);
s = s(:);
C = (s.*(log(s)/log(3)))\runtime;
h2 = plot(3.^(begin:last), C*s.*log(s)/log(3), '--', Color='r', LineWidth=1);
% h2 = plot(xx,yy,LineStyle='--',Color='r',LineWidth=0.8);
legend([h1,h2],'runtime',['y=',num2str(C),'*n*log_3n'],Location='southeast')



%%
function y = radix_3_fft(x)
    x = x(:);
    n = length(x);
    omega = exp(-2*pi*1i/n);        
    t1 = exp(-2*pi*1i/3);
    t2 = exp(-4*pi*1i/3);
    if rem(n,3) == 0
        k1 = (0:n/3-1)';
        k2 = 2*k1;
        w1 = omega.^k1;
        w2 = omega.^k2;
        u0 = radix_3_fft(x(1:3:n-2));
        u1 = w1.*radix_3_fft(x(2:3:n-1));
        u2 = w2.*radix_3_fft(x(3:3:n));
        y = [u0 + u1 + u2; 
            u0 + t1*u1 + t2*u2; 
            u0 + t2*u1 + t1*u2];
    else
        j = 0:n-1;
        k = j';
        F = omega.^(k*j);
        y = F*x;
    end
end
