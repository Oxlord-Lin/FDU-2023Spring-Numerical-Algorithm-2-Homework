clear
clc
%% 毕达哥拉斯
offset = 0;
x1 = 0; y1 = 0;
vx1 = 0; vy1 = 0;
x2 = 0; y2 = 4;
vx2 = 0; vy2 = offset;
x3 = 3; y3 = 0;
vx3 = 0; vy3 = 0;
planet1 = [x1,   vx1,   y1,  vy1];
planet2 = [x2,   vx2,   y2,  vy2];
planet3 = [x3,   vx3,   y3,  vy3];
t = [0,60];
ic = [planet1,planet2,planet3];
orbit(t,ic,t(2)*100);
pause(2)
tit = ['毕达哥拉斯问题，扰动\Delta v_{y_{2}}=',num2str(offset),'，时间范围：0~',num2str(t(2))];
title(tit)

% saveas(gcf,[tit,'.png'])
%%
function m1 = GetGlobal_m1()
% 用于获取第一个天体的质量
    m1 = 5;
end
function m2 = GetGlobal_m2()
% 用于获取第2个天体的质量
    m2 = 3;
end
function m3 = GetGlobal_m3()
% 用于获取第3个天体的质量
    m3 = 4;
end
function g = GetGlobal_g()
% 用于获取万有引力常数，这里归一化为1
    g = 1;
end


function z=orbit(inter,ic,n)
% input如下
% 时间区间 inter
% 初始条件 ic=[x1,vx1,y1,vy1,x2,vx2,y2,vy2,x3,vx3,y3,vy3]
% 分别表示3个天体的初始x坐标，初始x方向速度，初始y坐标，初始y方向速度
% 步数n
% output：绘制出天体的运行轨迹，并且返回最终的状态z
% 使用4阶龙格-库塔方法求解
    h = (inter(2)-inter(1))/n; % 步长
    y(1,:) = ic; %初始状态
    t(1) = inter(1);
%     set(gca,'XLim',[-3,3],'YLim',[-3,3],'XTick',[-3,0,3],'YTick',[-3,0,3],'SortMethod','childorder','Visible','on');
%     cla;
%     head1 = line('color','r','Marker','.','markersize',30,'xdata',[],'ydata',[]);
%     tail1 = line('color','r','LineStyle','-','xdata',[],'ydata',[]);
%     head2 = line('color','b','Marker','.','markersize',30,'xdata',[],'ydata',[]);
%     tail2 = line('color','b','LineStyle','-','xdata',[],'ydata',[]);
%     head3 = line('color','c','Marker','.','markersize',30,'xdata',[],'ydata',[]);
%     tail3 = line('color','c','LineStyle','-','xdata',[],'ydata',[]);
    
    % 4阶龙格-库塔法
    figure
    p = 10; % 每走10步就画图，这样节省空间，不用全部存储
    for k = 1:n/p
        for i = 1:p
            t(i+1) = t(i) + h;
            y(i+1,:) = RK4step(t(i),y(i,:),h);
        end     

        plot(y(:,1),y(:,3),'r:',LineWidth=0.75); 
        hold on
        plot(y(:,5),y(:,7),'b--',LineWidth=0.55);
        plot(y(:,9),y(:,11),'g-',LineWidth=0.55);
        y(1,:) = y(p+1,:); t(1) = t(p+1);
        t(end)
        pause(0.01)
    end
    
    head1 = scatter(y(1,1),y(1,3),30,'red','filled');
    head2 = scatter(y(1,5),y(1,7),30,'blue','filled');
    head3 = scatter(y(1,9),y(1,11),30,'g','filled');
    legend([head1,head2,head3],'天体1','天体2','天体3',Location='best')
    axis equal
    t = ['平面上的三体问题，','时间范围：0~',num2str(inter(2))];
    title(t)
%     saveas(gca,[t,'.png'])
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
    m3 = GetGlobal_m3();
    g = GetGlobal_g();
    mg1 = m1*g; mg2 = m2*g; mg3 = m3*g;
    x1 = Y(1); y1 = Y(3); 
    x2 = Y(5); y2 = Y(7);
    x3 = Y(9); y3 = Y(11);
    vx1 = Y(2); vy1 = Y(4); 
    vx2 = Y(6); vy2 = Y(8);
    vx3 = Y(10); vy3 = Y(12);
    dist12 = sqrt((x2-x1)^2 + (y2-y1)^2);
    dist13 = sqrt((x3-x1)^2 + (y3-y1)^2);
    dist23 = sqrt((x3-x2)^2 + (y3-y2)^2);

    z = zeros(1,12);
    z(1) = vx1;
    z(2) = (mg2*(x2-x1))/dist12^3 + (mg3*(x3-x1))/dist13^3;
    z(3) = vy1;
    z(4) = (mg2*(y2-y1))/dist12^3 + (mg3*(y3-y1))/dist13^3;

    z(5) = vx2;
    z(6) = (mg1*(x1-x2))/dist12^3 + (mg3*(x3-x2))/dist23^3;
    z(7) = vy2;
    z(8) = (mg1*(y1-y2))/dist12^3 + (mg3*(y3-y2))/dist23^3;

    z(9) = vx3;
    z(10)= (mg1*(x1-x3))/dist13^3 + (mg2*(x2-x3))/dist23^3;
    z(11) = vy3;
    z(12)= (mg1*(y1-y3))/dist13^3 + (mg2*(y2-y3))/dist23^3;
end