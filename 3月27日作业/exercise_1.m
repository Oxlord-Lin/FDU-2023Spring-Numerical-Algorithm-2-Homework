data=[-1.0000 -1.0000 1.6389
-1.0000 1.0000 0.5403
1.0000 -1.0000 -0.9900
1.0000 1.0000 0.1086
-0.7313 0.6949 0.9573
0.5275 -0.4899 0.8270
-0.0091 -0.1010 1.6936
0.3031 0.5774 1.3670];
x = data(:,1);
y = data(:,2);
z = data(:,3);
% Shepard's method
n = 100;
X = linspace(-1,1,n);
Y = X;
Z = zeros(n);
for i = 1:n
   for j = 1:n
       Z(i,j) = myfunc1(x,y,z,X(i),Y(j));
   end
end
% figure
mesh(X,Y,Z','FaceAlpha',0.25) % 这个转置至关重要 
hold on
scatter3(x,y,z,'red','filled','o');
xlabel('X')
ylabel('Y')
zlabel('Z')
title("2D interpolation with the Shepard's method")

%% piecewise interpolation based on Delaunay triangulization
rng default;
TR = delaunayTriangulation(x,y);
for i = 1:n
    for j = 1:n 
        P = [X(i),Y(j)];
        [ID,B] = pointLocation(TR,P); % B是相对于该点所属三角形的该点的重心坐标
        index = TR.ConnectivityList(ID,:); % index是包含该店的三角形的三个点的序号
        Z(i,j) = B(1)*z(index(1)) + B(2)*z(index(2)) + B(3)*z(index(3)); % z是input的第三列
    end
end
figure
mesh(X,Y,Z','FaceAlpha',0.25)
hold on
scatter3(x,y,z,'red','filled','o');
xlabel('X')
ylabel('Y')
zlabel('Z')
title("2D interpolation with the piecewise interpolation based on Delaunay Triangulization")


function Z = myfunc1(x,y,z,s,t)
    point = [s,t]';
    n = length(x);
    denom = 0; % 分母
    flg = 0;
    for i = 1:n
        denom = denom + (1/norm([x(i),y(i)]'-point));
        if isinf(denom)
            flg = 1; % 说明这里是插值点
            k = i;
            break
        end
    end
    Z = 0;
    for j = 1:n
        Z = Z + z(j)*((1/norm([x(j),y(j)]'-point))/denom);
    end
    if flg == 1
        Z = z(k);
    end
end