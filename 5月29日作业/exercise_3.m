exact_result = 0.5*(-2+exp(1)-sin(1)+cos(1));
error = [];
I_store = [];
for n = 1:50 % 分成n段
    y = linspace(1,0,n+1);
    x = linspace(0,1,n+1);
    I = 0;
    for i = 2:length(x)
        for j = 1:i-2
            I = I + tri_integration([x(j),y(i-1)],[x(j),y(i)],[x(j+1),y(i)]); % 下三角
            I = I + tri_integration([x(j+1),y(i)],[x(j+1),y(i-1)],[x(j),y(i-1)]); % 上三角
        end
        j = i-1;
        I = I + tri_integration([x(j),y(i-1)],[x(j),y(i)],[x(j+1),y(i)]); % 最右边的只有下三角
    end
    I_store = [I_store,I];
    error = [error,abs(I-exact_result)];
end
plot(1:n,I_store,'b-o',LineWidth=1);
xlabel('剖分段数n');
ylabel('数值积分结果');
title('用三角剖分进行数值积分的结果随剖分段数的变化')
saveas(gcf,'exercise-3-integration-result.png')
figure
semilogy(1:n,error,'r-.o',Linewidth=1);
xlabel('剖分段数n');
ylabel('与精确值的误差（取对数）');
title('数值积分的误差随三角剖分的大小的变化关系');
saveas(gcf,'exercise-3-integration-error.png')

function I = tri_integration(A,B,C)
    % xa = A(1); 
    ya = A(2);
    xb = B(1); 
    yb = B(2);
    xc = C(1); 
    % yc = C(2);
    S = 0.5*(ya-yb)^2;
    x1 = xb + 1/6 * (xc-xb);
    x2 = xb + 2/3 * (xc-xb);
    y1 = yb + 1/6 * (ya-yb);
    y2 = yb + 2/3 * (ya-yb);
    I = S* (1/3) * (fun(x1,y2) + fun(x2,y1) + fun(x1,y1));
end

function f = fun(x,y)
    f = exp(x)*sin(y);
end