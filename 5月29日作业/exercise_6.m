g1 = @(x) x + 1;
g2 = g1;
g3 = @(y) y^2 - 1;
g4 = @(y) y^2 + 1;

FDM_poisson(-1,1,-1,1,5,5,g1,g2,g3,g4);
% saveas(gcf,'exercise-6-5by5.png')

figure
FDM_poisson(-1,1,-1,1,10,10,g1,g2,g3,g4);
% saveas(gcf,'exercise-6-10by10.png')
saveas(gcf,'exercise-6-10by10-temp.png')

figure
FDM_poisson(-1,1,-1,1,20,20,g1,g2,g3,g4);
% saveas(gcf,'exercise-6-20by20.png')


% figure
% FDM_poisson(-1,1,-1,1,100,100,g1,g2,g3,g4);
% saveas(gcf,'exercise-6-100by100.png')

function w = FDM_poisson(xl,xr,yb,yt,M,N,g1,g2,g3,g4)
% input:矩形区域[xl,xr]*[yb,yt]，有M*N块，即M+1和N+1个均匀点
% 以及【下、上、左、右】四个边界的迪利克雷条件g1,g2,g3,g4(函数句柄)
% outpu:包含解的矩阵w，并作图
    m = M + 1;  n = N + 1;  mn = m*n;
    h = (xr - xl)/M;     h2 = h^2;
    k = (yt - yb)/N;     k2 = k^2;
    x = xl + (0:M)*h; % 设置网格值
    y = yb + (0:N)*k;
    A = zeros(mn); b = zeros(mn,1);
    % 中间部分
    P=[];
    for i = 2:m-1 % 中间部分
        for j = 2:n-1
            p = i + (j-1)*m;
            P = [P,p];
            A(p,i-1+(j-1)*m) = 1/h2;    A(p,i+1+(j-1)*m) = 1/h2; % 左右
            A(p,i+(j-1)*m) = -2/h2 - 2/k2;
            A(p,i+(j-2)*m) = 1/k2;      A(p,i+j*m) = 1/k2; % 下上
            % b(p) = f(x(i),y(j));
        end
    end
    % 上下边界
    for i = 1:m
        j = 1; A(i+(j-1)*m, i+(j-1)*m) = 1; b(i+(j-1)*m) = g1(x(i)); %下边界
        j = n; A(i+(j-1)*m, i+(j-1)*m) = 1; b(i+(j-1)*m) = g2(x(i)); %上边界
    end
    % 左右边界
    for j = 2:n-1
        i = 1; A(i+(j-1)*m, i+(j-1)*m) = 1; b(i+(j-1)*m) = g3(y(j)); %左边界
        i = m; A(i+(j-1)*m, i+(j-1)*m) = 1; b(i+(j-1)*m) = g4(y(j)); %右边界
    end
    v = A\b;
%     v = jacobi(A,b); % 稀疏矩阵，需要特殊处理，否则解不出来
    w = reshape(v(1:mn),m,n);
%     surf(x,y,w',LineStyle="none");
    surf(x,y,w',LineStyle="-");
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['有限差分法解拉普拉斯方程，','网格大小：',num2str(M),'\times',num2str(N)]);
end