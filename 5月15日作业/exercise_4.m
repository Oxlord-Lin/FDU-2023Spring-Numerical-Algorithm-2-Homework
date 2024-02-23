clear;
clc;
% 生成不同的步长
t = 4; % 一共要尝试几种步长
step_size = ones(1,t);
s = 0.5;
k = 0:t-1;
s = s.^k;
step_size = s.*step_size;

%% 不同步长算出的一阶导数 
for i = 1:length(step_size)
    % 生成数据点
    step = step_size(i);
    x = 1:step:4;
    y = func1(x);
    % 样条插值
    pp = spline(x,y); % 返回一个分段多项式结构体
    C = pp.coefs;
    [m,n] = size(C);
    % 计算一阶导数并绘图
    for j = 1:m
        xx = linspace(1+step*(j-1), 1+step*j, 500);
%         real_first_dev = exp(xx) + 1./xx;
%         real_second_dev = exp(xx) - 1./(xx.^2);
        x1 = x(j);
        coef = C(j,:);
        a = coef(1);
        b = coef(2);
        c = coef(3);
        d = coef(4);
        estimate_first_dev = 3*a*(xx-x1).^2 + 2*b*(xx-x1) + c;
%         estimate_second_dev = 6*a*(xx-x1) + 2*b;
        if i == 1
            h1 = plot(xx,estimate_first_dev,'r-',LineWidth=1); hold on
        elseif i == 2
            h2 = plot(xx,estimate_first_dev,'k-.',LineWidth=1.75); hold on
        elseif i == 3
            h3 = plot(xx,estimate_first_dev,'b:',LineWidth=1.5); hold on
        elseif i == 4
            h4 = plot(xx,estimate_first_dev,'m--',LineWidth=1); hold on
        end
    end
end
legend([h1,h2,h3,h4], ...
    ['step = ',num2str(step_size(1))], ...
    ['step = ',num2str(step_size(2))], ...
    ['step = ',num2str(step_size(3))], ...
    ['step = ',num2str(step_size(4))],Location='best');
title('插值样条函数的一阶导')
saveas(gcf,'exercise-4 一阶导.png')

%% 不同步长算出的二阶导数 
figure
for i = 1:length(step_size)
    % 生成数据点
    step = step_size(i);
    x = 1:step:4;
    y = func1(x);
    % 样条插值
    pp = spline(x,y); % 返回一个分段多项式结构体
    C = pp.coefs;
    [m,n] = size(C);
    % 计算二阶导数并绘图
    for j = 1:m
        xx = linspace(1+step*(j-1), 1+step*j, 500);
%         real_first_dev = exp(xx) + 1./xx;
%         real_second_dev = exp(xx) - 1./(xx.^2);
        x1 = x(j);
        coef = C(j,:);
        a = coef(1);
        b = coef(2);
        c = coef(3);
        d = coef(4);
%         estimate_first_dev = 3*a*(xx-x1).^2 + 2*b*(xx-x1) + c;
        estimate_second_dev = 6*a*(xx-x1) + 2*b;
        if i == 1
            h1 = plot(xx,estimate_second_dev,'r-',LineWidth=1); hold on
        elseif i == 2
            h2 = plot(xx,estimate_second_dev,'k-.',LineWidth=1.75); hold on
        elseif i == 3
            h3 = plot(xx,estimate_second_dev,'b:',LineWidth=1.5); hold on
        elseif i == 4
            h4 = plot(xx,estimate_second_dev,'m--',LineWidth=1); hold on
        end
    end
end
legend([h1,h2,h3,h4], ...
    ['step = ',num2str(step_size(1))], ...
    ['step = ',num2str(step_size(2))], ...
    ['step = ',num2str(step_size(3))], ...
    ['step = ',num2str(step_size(4))],location='best');
title('插值样条函数的二阶导')
saveas(gcf,'exercise-4 二阶导.png')


%% 不同步长算出的一阶导数的误差 
figure
for i = 1:length(step_size)
    % 生成数据点
    step = step_size(i);
    x = 1:step:4;
    y = func1(x);
    % 样条插值
    pp = spline(x,y); % 返回一个分段多项式结构体
    C = pp.coefs;
    [m,n] = size(C);
    % 计算一阶导数的误差并绘图
    for j = 1:m
        xx = linspace(1+step*(j-1), 1+step*j, 500);
        real_first_dev = exp(xx) + 1./xx;
%         real_second_dev = exp(xx) - 1./(xx.^2);
        x1 = x(j);
        coef = C(j,:);
        a = coef(1);
        b = coef(2);
        c = coef(3);
        d = coef(4);
        estimate_first_dev = 3*a*(xx-x1).^2 + 2*b*(xx-x1) + c;
%         estimate_second_dev = 6*a*(xx-x1) + 2*b;

%         if i == 1
%             h1 = semilogy(xx,abs(estimate_first_dev - real_first_dev),'r-'); hold on
%         elseif i == 2
%             h2 = semilogy(xx,abs(estimate_first_dev - real_first_dev),'k-.'); hold on
%         elseif i == 3
%             h3 = semilogy(xx,abs(estimate_first_dev - real_first_dev),'b:'); hold on
%         elseif i == 4
%             h4 = semilogy(xx,abs(estimate_first_dev - real_first_dev),'m--'); hold on
%         end

        if i == 1
            h1 = plot(xx,log2(abs(estimate_first_dev - real_first_dev)),'r-'); hold on
        elseif i == 2
            h2 = plot(xx,log2(abs(estimate_first_dev - real_first_dev)),'k-.'); hold on
        elseif i == 3
            h3 = plot(xx,log2(abs(estimate_first_dev - real_first_dev)),'b:'); hold on
        elseif i == 4
            h4 = plot(xx,log2(abs(estimate_first_dev - real_first_dev)),'m--'); hold on
        end
    end

end
legend([h1,h2,h3,h4], ...
    ['step = ',num2str(step_size(1))], ...
    ['step = ',num2str(step_size(2))], ...
    ['step = ',num2str(step_size(3))], ...
    ['step = ',num2str(step_size(4))],location='best');
title('插值样条函数的一阶导误差（取对数）')
ylim([-21,3])
yticks([-21,-18,-15,-12,-9,-6,-3,0,3])
yticklabels({'2^{-21}','2^{-18}','2^{-15}','2^{-12}','2^{-9}','2^{-6}','2^{-3}','2^{0}','2^{3}'})
saveas(gcf,'exercise-4 一阶导的误差.png')

%% 不同步长算出的二阶导数的误差 
figure()
for i = 1:length(step_size)
    % 生成数据点
    step = step_size(i);
    x = 1:step:4;
    y = func1(x);
    % 样条插值
    pp = spline(x,y); % 返回一个分段多项式结构体
    C = pp.coefs;
    [m,n] = size(C);
    % 计算二阶导数的误差并绘图
    for j = 1:m
        xx = linspace(1+step*(j-1), 1+step*j, 500);
%         real_first_dev = exp(xx) + 1./xx;
        real_second_dev = exp(xx) - 1./(xx.^2);
        x1 = x(j);
        coef = C(j,:);
        a = coef(1);
        b = coef(2);
        c = coef(3);
        d = coef(4);
%         estimate_first_dev = 3*a*(xx-x1).^2 + 2*b*(xx-x1) + c;
        estimate_second_dev = 6*a*(xx-x1) + 2*b;

%         if i == 1
%             h1 = semilogy(xx,abs(estimate_second_dev - real_second_dev),'r-'); hold on
%         elseif i == 2
%             h2 = semilogy(xx,abs(estimate_second_dev - real_second_dev),'k-.'); hold on
%         elseif i == 3
%             h3 = semilogy(xx,abs(estimate_second_dev - real_second_dev),'b:'); hold on
%         elseif i == 4
%             h4 = semilogy(xx,abs(estimate_second_dev - real_second_dev),'m--'); hold on
%         end

         if i == 1
            h1 = plot(xx,log2(abs(estimate_second_dev - real_second_dev)),'r-'); hold on
        elseif i == 2
            h2 = plot(xx,log2(abs(estimate_second_dev - real_second_dev)),'k-.'); hold on
        elseif i == 3
            h3 = plot(xx,log2(abs(estimate_second_dev - real_second_dev)),'b:'); hold on
        elseif i == 4
            h4 = plot(xx,log2(abs(estimate_second_dev - real_second_dev)),'m--'); hold on
        end
    end
end
legend([h1,h2,h3,h4], ...
    ['step = ',num2str(step_size(1))], ...
    ['step = ',num2str(step_size(2))], ...
    ['step = ',num2str(step_size(3))], ...
    ['step = ',num2str(step_size(4))],location='best');
title('插值样条函数的二阶导误差（取对数）')
ylim([-16,5])
yticks([-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0,2,4])
yticklabels({'2^{-20}','2^{-18}','2^{-16}','2^{-14}','2^{-12}','2^{-10}','2^{-8}','2^{-6}','2^{-4}','2^{-2}','2^{0}','2^{2}','2^{4}'})
saveas(gcf,'exercise-4 二阶导的误差.png')


function y = func1(x)
    n = length(x);
    y = zeros(1,n);
    for i = 1:n
        y(i) = exp(x(i)) + log(x(i));
    end
end