T = 2;
for d = 6:2:12
    x = linspace(-1,1,d+1)';
    f = myfunc(x(1:end-1));
    coef = 1/d*fft(f); % coef是系数
%     a = real(coef);
%     b = imag(coef);
    % 绘图
    t = 10000;
    % 这个方法画出来的图不正确
    xx = linspace(-1,1,t);
    S = zeros(t,1);
    for j = 1:t
      S(j) = 0;
      for k = 0:d-1
        S(j) = S(j) + coef(k+1)*exp(-1i*2*pi*k*xx(j)/T);
      end
    end
    yy = real(S);
    figure
    h1 = plot(xx,yy); 
    hold on
    h2 = plot(x(1:end-1),f,'*');
    legend([h1,h2],[num2str(d),'阶三角插值多项式'],'选取的插值点')
    
end


function f = myfunc(x)
    f = zeros(length(x),1);
    for i = 1:length(x)
        if x(i) == -1 || x(i) == 0 ||x(i) == 1
            f(i) = 0;
        elseif x(i)<0
                f(i) = -1;
        else 
            f(i) = 1;
        end
    end
end