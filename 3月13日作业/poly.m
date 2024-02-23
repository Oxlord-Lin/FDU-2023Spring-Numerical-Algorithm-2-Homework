for k = 3:13
    node = linspace(-1,1,k);
    t = linspace(-1,1,1000);
    f = ones(1,1000);
    for i = 1:1000
        v = t(i);
        for j = 1:k
            f(i) = f(i)*abs(v-node(j));
        end
    end
    plot(t,f,'DisplayName',['节点数=',num2str(k)],LineWidth=1.5);
    hold on
end
legend