nodes = linspace(0.2,0.8,4); % 初始猜测的交错点
f = log(nodes+1)'; % 猜测节点处的函数值
k = ones(4,1); 
epsilon = 1e-16; % 精度
for i = 1:100 % 迭代保护
    A = myfunc1(nodes);
    k = A\f;  % k = [a,b,c, eta]' a,b,c 分别是二次多项式的系数，eta是节点上的误差 f(x)-p(x)
    a = k(1);
    b = k(2);
    c = k(3);
    eta = k(4);
    xx = linspace(0,1,1001);
    f = @(x) log(1+x) - a*x.^2 - b*x - c;
    err = f(xx);
    % 绘制误差曲线
    if i == 1
        ax = nexttile;
        plot(ax,xx,err);
        hold on
        scatter(ax,nodes,f(nodes),'o','filled')
        [~,~] = title('the initial alternating set and the error curve',['a=',num2str(-a),',b=',num2str(-b),',c=',num2str(-c)]);
%         legend('error curve','alternating set')
    end
    if i == 2 || i == 3
        ax = nexttile;
        plot(ax,xx,err);
        hold on
        scatter(ax,nodes,f(nodes),'o','filled')
        [~,~] = title(['the error curve of loop=',num2str(i-1)],['a=',num2str(-a),',b=',num2str(-b),',c=',num2str(-c)]);
%         legend('error curve','alternating set')
    end
    [max_abs_err, max_index] = max(abs(err)); % max_abs_err 是 f(x)-p(x) 的无穷范数

    if abs(max_abs_err - abs(eta)) < epsilon
        break
    end

    % 使用【单一交换法】进行结点更新
    x_new = xx(max_index); % 将要插入的新节点
    max_err = err(max_index); % 新节点处误差
    new_nodes = update_nodes(x_new,max_err,nodes,eta); % 插入新节点
    nodes = new_nodes;
    f = log(nodes+1)';
end
% 绘图，最终结果
ax = nexttile;
plot(ax,xx,err);
hold on
scatter(ax,nodes,f(nodes),'o','filled')
[~,~] = title(['the final error curve of loop=',num2str(i-1)],['a=',num2str(-a),',b=',num2str(-b),',c=',num2str(-c)]);
legend('error curve','alternating set',Location='bestoutside')

function A = myfunc1(nodes) % 某一组交错点所对应的矩阵
    x1 = nodes(1);
    x2 = nodes(2);
    x3 = nodes(3);
    x4 = nodes(4);
    A = [x1^2, x1, 1, 1;
         x2^2, x2, 1, -1;
         x3^2, x3, 1, 1;
         x4^2, x4, 1, -1];
end

function new_nodes = update_nodes(x_new,max_err,nodes,eta)
    flg = 0;
    if x_new <= nodes(1)
        if max_err * eta > 0  % 如果与x1处误差同号
            nodes(1) = x_new; % 用x_new替换x1
            flg = 1;
        else                %如果与x1处误差异号
            nodes(2:4) = nodes(1:3); %把原来最后一个节点扔掉
            nodes(1) = x_new;
            flg = 1;
        end    
    end


    for j = 1:3
        if flg == 1
            break 
        end
        if nodes(j) < x_new && x_new <= nodes(j+1)
            if max_err * eta * (-1)^(j-1) > 0 % 如果与x(j)处误差同号
                nodes(j) = x_new; % 用x_new替换x(j)
                flg = 1;
                break % 如果与x(j+1)处误差同号
            else % 如果与x(j+1)处误差同号
                nodes(j+1) = x_new; % 用x_new替换x(j+1)
                flg = 1;
                break
            end
        end
    end

    if x_new > nodes(4)
        if max_err *(-1) *eta > 0 % 如果与x4处误差同号
            nodes(4) = x_new; % 用x_new替换x4
        else                    % 如果与x4处误差异号
            nodes(1:3) = nodes(2:4); % 把原来第一个节点扔掉
            nodes(4) = x_new;
        end    
    end
    new_nodes = nodes;
end