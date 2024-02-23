clear
clc
t_final = 20000;
ic = [2, 0.2, 2, -0.2,  0, -0.01, 0, 0.01];

%% 
orbit([0,t_final],ic,t_final*100,1);
%%
orbit([0,t_final],ic,t_final*100,2);
%%
orbit([0,t_final],ic,t_final*100,3);
% 1：显式梯形，2：Störmer-Verlet，3：RK4

function g = GetGlobal_g()
% 用于获取万有引力常数，这里归一化为1
    g = 1;
end
function m1 = GetGlobal_m1()
% 用于获取第一个天体的质量
    m1 = 0.3;
end
function m2 = GetGlobal_m2()
% 用于获取第2个天体的质量
    m2 = 0.03;
end
function z=orbit(inter,ic,n,method)
% input如下
% 时间区间 inter
% 初始条件 ic=[x1,vx1,y1,vy1,x2,vx2,y2,vy2]
% 分别表示两个天体的初始x坐标，初始x方向速度，初始y坐标，初始y方向速度
% 步数n
% 调用单步方法,用method指定
% output：绘制出两个天体的运行轨迹，并且返回最终的坐标z
% 1：显式梯形，2：Störmer-Verlet，3：4阶龙格-库塔
    h = (inter(2)-inter(1))/n; % 步长
    y = ic; %初始状态
    E0 = total_energy(y(1,:)); % 初始总能量
    energy_diff_ex_trap=[E0];
    energy_diff_SV = energy_diff_ex_trap;
    energy_diff_RK4 = energy_diff_ex_trap;
    ts=[0];
    t(1) = inter(1);
%     set(gca,'XLim',[-3,3],'YLim',[-3,3],'XTick',[-3,0,3],'YTick',[-3,0,3],'SortMethod','childorder','Visible','on');
%     cla;
%     head1 = line('color','r','Marker','.','markersize',30,'xdata',[],'ydata',[]);
%     tail1 = line('color','r','LineStyle','-','xdata',[],'ydata',[]);
%     head2 = line('color','b','Marker','.','markersize',30,'xdata',[],'ydata',[]);
%     tail2 = line('color','b','LineStyle','--','xdata',[],'ydata',[]);

    % 显式梯形法
    if method == 1
        for k = 1:n
            t = t + h;
            y = trapstep(t,y,h);
            if mod(k,4000) == 3999
%                 energy_diff_ex_trap = [energy_diff_ex_trap, -E0 + total_energy(y(k+1,:)) ];
                energy_diff_ex_trap = [energy_diff_ex_trap, total_energy(y) ];
                ts = [ts,t];
            end
        end
        z = y;
        z = z(1:2:end);
    % Störmer-Verlet法
    elseif method == 2
        v = ic(2:2:end);
        r = ic(1:2:end);
        v = SV_v_step(t(1),r,v,0.5*h); % 第一次，速度走"半步"
        for k = 1:n
            t = t + h;
            r = r + h*v;
            v = SV_v_step(t,r,v,h);
            if mod(k,20000) == 19999
                Y_temp = zeros(1,8);
                Y_temp(1:2:7) = r;
                Y_temp(2:2:8) = v; % 只能近似一下，毕竟得不得到同一个时间点的速度和位置
%                 energy_diff_SV = [energy_diff_SV, -E0 + total_energy(Y_temp) ];
                energy_diff_SV = [energy_diff_SV,  total_energy(Y_temp) ];
                ts = [ts,t];
            end
        end
        z = y;
        z = z(1:2:end);
    % 4阶龙格-库塔法
    elseif method == 3
        for k = 1:n
            t = t + h;
            y = RK4step(t,y,h);
            if mod(k,1000) == 999
%                 energy_diff_RK4 = [energy_diff_RK4, -E0 + total_energy(y(k+1,:)) ];
                energy_diff_RK4 = [energy_diff_RK4, total_energy(y) ];
                ts = [ts,t];
            end
        end
        z = y;
        z = z(1:2:end);
    end
%     if method ~= 2
%         set(head1,'xdata',y(n,1),'ydata',y(n,3));
%         set(tail1,'xdata',y(1:n-1,1),'ydata',y(1:n-1,3));
%         set(head2,'xdata',y(n,5),'ydata',y(n,7));
%         set(tail2,'xdata',y(1:n-1,5),'ydata',y(1:n-1,7));
%     elseif method == 2
%         set(head1,'xdata',r(n,1),'ydata',r(n,2));
%         set(tail1,'xdata',r(1:n-1,1),'ydata',r(1:n-1,2));
%         set(head2,'xdata',r(n,3),'ydata',r(n,4));
%         set(tail2,'xdata',r(1:n-1,3),'ydata',r(1:n-1,4));
%     end
%     drawnow;
%     legend([head1,head2],'天体1','天体2',Location='best')
%     axis equal
%     if method == 1
%         t = ['平面上的二体问题,','显式梯形法,','步长=',num2str(h),',时间范围：0~',num2str(inter(2))];
%         title(t)
%     elseif method == 2
%         t = ['平面上的二体问题,',' Störmer-Verlet法,','步长=',num2str(h),',时间范围：0~',num2str(inter(2))];
%         title(t)
%     elseif method == 3
%         t = ['平面上的二体问题,',' 4阶龙格-库塔法,','步长=',num2str(h),',时间范围：0~',num2str(inter(2))];
%         title(t)
%     end
%     saveas(gca,[t,'.png'])
    
    figure
    if method == 1
        semilogy(ts,energy_diff_ex_trap,'m-.',LineWidth=0.5);
        xlabel('时间')
        ylabel('机械能')
        t = ['显式梯形法的机械能变化,','时间范围：0~',num2str(inter(2))];
    elseif method == 2
        semilogy(ts,energy_diff_SV,'b-',LineWidth=0.5)
        xlabel('时间')
        ylabel('机械能（近似值）')
        t = ['Störmer-Verlet法的机械能（近似值）变化,','时间范围：0~',num2str(inter(2))];
    elseif method == 3
        semilogy(ts,energy_diff_RK4,'k--',LineWidth=1)
        xlabel('时间')
        ylabel('当前总能量')
        t = ['4阶龙格-库塔法的总能量变化,','时间范围：0~',num2str(inter(2))];
    end
    title(t);
%     ylim(1e-3*[8.82,8.822])
    saveas(gcf,[t,'.png'])
end

function y = trapstep(t,x,h)
% 显示梯形法
% input:x是当前状态，t是当前时间，h是步长
% output：y是下一刻的状态
    z1 = ydot(t,x);
    g = x+h*z1;
    z2 = ydot(t+h,g);
    y = x + h*(z1 + z2)/2;
end

function v_next = SV_v_step(t,r,v,h)
% Störmer-Verlet法，本函数用于更新v
% input:r是k步位置，v是(k-1/2)步的速度，t是当前时间，h是步长
% output：v_next是(k+1/2)步的速度
    m1 = GetGlobal_m1();
    m2 = GetGlobal_m2();
    g = GetGlobal_g();
    mg1 = m1*g; mg2 = m2*g;
    x1 = r(1); y1 = r(2); x2 = r(3); y2 = r(4);
    vx1 = v(1); vy1 = v(2); vx2 = v(3); vy2 = v(4);
    dist = sqrt((x2-x1)^2 + (y2-y1)^2); % 两个星体之间的距离
    dist3 = dist^3;
    v_next = zeros(1,4);
    v_next(1) = vx1 + h * (-mg2*(x1-x2)) / dist3;
    v_next(2) = vy1 + h * (-mg2*(y1-y2)) / dist3;
    v_next(3) = vx2 + h * (-mg1*(x2-x1)) / dist3;
    v_next(4) = vy2 + h * (-mg1*(y2-y1)) / dist3;
end

function y = RK4step(t,x,h)
% 4阶龙格-库塔法
% input:x是当前状态，t是当前时间，h是步长
% output：y是下一刻的状态
    s1 = ydot(t,x);
    s2 = ydot(t+h/2,x+h/2*s1);
    s3 = ydot(t+h/2,x+h/2*s2);
    s4 = ydot(t+h,x+h*s3);
    y = x + h/6 * (s1 + 2*s2 + 2*s3 + s4);
end

function z = ydot(t,Y)
% input:x是当前状态，t是当前时间，h是步长
% output：z=ydot，是x处的导数，z=dY/dt
    m1 = GetGlobal_m1();
    m2 = GetGlobal_m2();
    g = GetGlobal_g();
    mg1 = m1*g; mg2 = m2*g;
    x1 = Y(1); y1 = Y(3); x2 = Y(5); y2 = Y(7);
    vx1 = Y(2); vy1 = Y(4); vx2 = Y(6); vy2 = Y(8);
    dist = sqrt((x2-x1)^2 + (y2-y1)^2);
    dist3 = dist^3;

    z = zeros(1,8);
    z(1) = vx1;
    z(2) = (mg2*(x2-x1))/dist3;
    z(3) = vy1;
    z(4) = (mg2*(y2-y1))/dist3;
    z(5) = vx2;
    z(6) = (mg1*(x1-x2))/dist3;
    z(7) = vy2;
    z(8) = (mg1*(y1-y2))/dist3;
end

function E = total_energy(Y)
% 求出Y状态下的系统总机械能
    m1 = GetGlobal_m1;
    m2 = GetGlobal_m2;
    g = GetGlobal_g;
    x1 = Y(1); y1 = Y(3); x2 = Y(5); y2 = Y(7);
    vx1 = Y(2); vy1 = Y(4); vx2 = Y(6); vy2 = Y(8);
    dist = sqrt((x2-x1)^2 + (y2-y1)^2);
    E = 0.5*m1*(vx1^2 + vy1^2) + 0.5*m2*(vx2^2 + vy2^2) - g*m1*m2/dist;
end