data = [1.00000 0.00000 -1.0000;
        0.80902 0.58779 -2.6807;
        0.30902 0.95106 5.6161;
        -0.30902 0.95106 5.6161;
        -0.80902 0.58779 -2.6807;
        -1.00000 0.00000 -1.0000;
        -0.80902 -0.58779 -2.6807;
        -0.30902 -0.95106 5.6161;
        0.30902 -0.95106 5.6161;
        0.80902 -0.58779 -2.6807];
a = data(:,1);
b = data(:,2);
i = sqrt(-1);
input = a + i*b;
output = data(:,3);
coef = ceof_Newton_basis(input,output);
% c = ones(10,1);
% for j = 1:10
%     for k = 1:10
%         if k == j 
%            continue
%         else
%             c(j) = c(j)*(input(k)-input(j));
%         end
%         c(j) = 1/c(j);
%     end
% end
%% 在单位圆上进行插值
t = linspace(0,2*pi,1000);
x = cos(t);
y = sin(t);
z = zeros(1000,1);
for k = 1:1000
    z(k) = real(interpolation(input,cos(t(k))+i*sin(t(k)),coef));
end
plot3(x,y,z,LineWidth=1.5)
hold on
plot3(x,y,zeros(1000),'--',Color='g')
xlabel('x轴，作为实部')
ylabel('y轴，作为虚部')
zlabel('z轴，只取实部')
scatter3(a,b,output,'filled','r')
title('在复数域上的插值多项式效果图（限制在单位圆上，取实部）')

%% 感觉不太成功
% X = [-20:20];
% Y = X;
% Z = zeros(1);
% for p = 1:41
%     for q = 1:41
%         Z(p,q) = real(interpolation(input,X(p)+i*Y(q),coef));
%     end
% end

step = .01;
s  = -1:step:1;
t = s;
[A,B] = meshgrid(s,t);
C = zeros(length(s),length(t));
for i = 1:length(s)
    for j = 1:length(t)
        x = A(i,j);
        y = B(i,j);
        if x^2 + y^2 >1;
            C(i,j) = 0;
            continue
        end
        C(i,j) = real(interpolation(input,x+1i*y,coef));
    end
end


figure()
colormap("default")
surf(A,B,C,'FaceAlpha',0.8);
view([20,20,120]);
shading interp


% surf(X,Y,Z,'FaceAlpha',0.5);
% hold on

colorbar;
xlabel('x')
ylabel('y')
zlabel('z')
% scatter3(a,b,output,'filled','r')
title('复平面上插值多项式示意图')

%%
% function f_z = myfunc(input,output,c,z) % 重心形式的Lagrange插值
%     L = zeros(10,1);
%     R_z = 1;
%     for j = 1:10
%         R_z = R_z*(z-input(j));
%     end
%     for k = 1:10
%         L(k) = R_z/(z-input(k))*c(k);
%     end
%     f_z = output'*L;
% end

function coef = ceof_Newton_basis(x,y)    
    y = y';
    x = x';
    n = length(x);
    divided_difference = zeros(n,n);
    divided_difference(1,:) = y; % 第一行作为0阶差商
    for i = 2:n
        for j = 1:n+1-i
            divided_difference(i,j) = (divided_difference(i-1,j+1) ...
                -divided_difference(i-1,j))/(x(j+i-1)-x(j));
        end
    end
    coef = divided_difference(:,1); % 牛顿插值的系数；
end

function f = interpolation(input,x,coef) % 计算插值多项式在该点的函数值
n = length(input);
temp = ones(n,1);
for l = 2:n
    temp(l) = temp(l-1)*(x-input(l-1));
end
f = coef'*temp;
end