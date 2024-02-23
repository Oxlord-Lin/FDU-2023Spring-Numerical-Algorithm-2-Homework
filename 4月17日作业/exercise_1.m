f = @(x) 1/sqrt((1+x)); % the original function
R11 = @(x) (1 + 0.5*x) / (1 + x); % R1/1 Pade approximant
R12 = @(x) (1 + 1.5*x) / (1 + 2*x + 0.5*x^2); % R1/2 Pade approximant
R21 = @(x) (1 - 0.75*x^2) / (1 + 0.5*x); % R2/1 Pade approximant

xx = linspace(0,2,2000);
yy11 = zeros(2000,1);
yy12 = yy11;
yy21 = yy11;

for i = 1:2000
    temp = xx(i);
    yy11(i) = R11(temp) - f(temp);
    yy12(i) = R12(temp) - f(temp);
    yy21(i) = R21(temp) - f(temp);
end

figure
plot(xx,yy11)
title('the approximation error of $R_{11}(x)=\frac{1+0.5x}{1+x}$','Interpreter','latex')
figure
plot(xx,yy12)
title('the approximation error of $R_{12}(x)=\frac{1+1.5x}{1+2x+0.5x^2}$','Interpreter','latex')
figure
plot(xx,yy21)
title('the approximation error of $R_{21}(x)=\frac{1-0.75x^2}{1+0.5x}$','Interpreter','latex')