x = [1,1,1,2,2,2,3,3];
Y = [3,3,3,1,1,1,2,2;
    3,3,NaN,1,1,NaN,2,0;
    3,NaN,NaN,1,NaN,NaN,0,0;
    NaN*ones(5,8)]; %NaN 表示尚待填入
coef = ceof_Newton_basis(x,Y); % coef中存储了Hermite插值的各阶系数
function coef = ceof_Newton_basis(x,Y)
n = length(x);
    for i = 2:n
        for j = 1:n+1-i
            if ~isnan(Y(i,j)) %已经提前填入
                Y(i,j) = Y(i,j);
            else
            Y(i,j) = (Y(i-1,j+1) - Y(i-1,j)) / (x(j+i-1) - x(j));
            end
        end
    end
    coef = Y(:,1); % 牛顿插值的系数；
end