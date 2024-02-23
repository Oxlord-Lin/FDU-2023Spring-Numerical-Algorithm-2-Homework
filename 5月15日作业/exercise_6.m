n = 3; % 初始的n
k = 4; % 外插的阶数为2*k
n_max = 2^(k-1) * n; % 最大的n
inter_results = Richardson_extrapolation(n,k);
plot(n*2.^(0:k-1),inter_results,'b-o');
xlabel('外插过程中用到的最大的n')
ylabel('用外插计算的\pi的近似值')
title('利用Richardson外插法高效计算pi——向刘徽与祖冲之的工作致敬')
format long
n*2.^(0:k-1)
inter_results
function inter_results = Richardson_extrapolation(n,k)
    Table = zeros(k,k);
    for i = 1:k
        Table(i,1) = myfunc(2^(i-1) * n);
    end
    for j = 2:k
        for i = 1:k-j+1
            Table(i,j) = ( 2^(2*j-2)*Table(i+1,j-1) - Table(i,j-1) ) / (2^(2*j-2) - 1);
        end
    end
    inter_results = Table(1,:);
end

function f = myfunc(n)
    f = n*sin(pi/n);
end