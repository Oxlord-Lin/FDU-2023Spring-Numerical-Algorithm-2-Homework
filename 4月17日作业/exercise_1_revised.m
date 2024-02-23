f = @(x) 1/sqrt((1+x)); % the original function
R11 = @(x) (1 + 0.25*x) / (1 + 0.75*x); % R1/1 Pade approximant
R02 = @(x) 1 / (1 + 0.5*x - 0.125*x^2); % R0/2 Pade approximant
R12 = @(x) (1 + 0.5*x) / (1 + x + 0.125*x^2); % R1/2 Pade approximant

x = linspace(0,2,2000);
Y = x;
y11 = zeros(2000,1);
y02 = yy11;
y12= yy11;

for i = 1:2000
    temp = x(i);
    Y(i) = f(temp);
    y11(i) = R11(temp);
    y02(i) = R02(temp);
    y12(i) = R12(temp);
end

figure
h1 = plot(x,y11,LineWidth=1); 
hold on
h2 = plot(x,Y,'--',LineWidth=1);
legend([h1,h2],'R_{1/1}(x)','1/sqrt(1+x)')
title('the Pade approximant $R_{1/1}(x)=\frac{1+0.25x}{1+0.75x}$','Interpreter','latex')
figure
h3 = plot(x,y02,LineWidth=1); hold on
h4 = plot(x,Y,'--',LineWidth=1);
legend([h3,h4],'R_{0/2}(x)','1/sqrt(1+x)')
title('the Pade approximant $R_{0/2}(x)=\frac{1}{1+0.5x-0.125x^2}$','Interpreter','latex')
figure
h5 = plot(x,y12); hold on
h6 = plot(x,Y,'--',LineWidth=1);
legend([h5,h6],'R_{1/2}(x)','1/sqrt(1+x)')
title('the Pade approximant $R_{1/2}(x)=\frac{1+0.5x}{1+x+0.125x^2}$','Interpreter','latex')