clear;clc;
T = 2;
for n = 6:6:24
    M = 2*n+1;
    x = linspace(0,2,M+1);
    f = myfunc(x(1:end-1));
    j1 = 0:M-1;
    j2 = j1 - n;
    coef = 1/M*fft(f);
    coef = fftshift(coef);

    t = 1000;
    xx = linspace(0,2,t);
    S = zeros(t,1);
    for j = 1:t
        S(j) = exp(1i*(xx(j)*j2)*2*pi/T)*coef;
    end
    yy = real(S);
    figure
    h1 = plot(xx,yy); 
    hold on
    h2 = plot(x(1:end-1),f,'*');
    legend([h1,h2],[num2str(n),'阶三角插值多项式'],'选取的插值点')
end


function f = myfunc(x)
    f = zeros(length(x),1);
    for i = 1:length(x)
        if x(i) == -1 || x(i) == 0 ||x(i) == 1 || x(i) == 2
            f(i) = 0;
        elseif x(i)<0 || x(i)>1
                f(i) = -1;
        else 
            f(i) = 1;
        end
    end
end