clear
clc
%% 限制性三体问题 漂亮的周期性轨道
% Broucke R 1  
% http://three-body.ipb.ac.rs/bsol.php?id=16
x1 = 0.8083106230; y1 = 0;
vx1 = 0; vy1 = 0.9901979166;
x2 = -0.4954148566; y2 = 0;
vx2 = 0; vy2 = -2.7171431768;
x3 = -0.3128957664; y3 = 0;
vx3 = 0; vy3 = 1.7269452602;

% % Henon 35 
% http://three-body.ipb.ac.rs/hsol.php?id=33
% x1= 0.0797756841 ;y1= 0;
% x2= 1.1140666180 ;y2= 0;
% x3= -1.1938423022 ;y3=-0;
% vx1= 0; vy1= 1.1801414216;
% vx2= 0; vy2= -0.2727205239;
% vx3= 0; vy3= -0.9074208977;

% % Oval, catface, and starship   需要高精度小步长
% http://three-body.ipb.ac.rs/sheen_sol.php?id=1 
% % 需要高精度，小步长
% x1= 0.536387073390 ;y1= 0.054088605008;
% x2= -0.252099126491 ;y2= 0.694527327749;
% x3= -0.275706601688 ;y3=-0.335933589318;
% vx1= -0.569379585581; vy1= 1.255291102531;
% vx2= 0.079644615252 ;vy2= -0.458625997341;
% vx3= 0.489734970329; vy3= -0.796665105189;

% GOGGLES  需要高精度小步长
% http://three-body.ipb.ac.rs/sol.php?id=6
% x1= -1 ;y1= 0;
% x2= 1 ;y2= 0;
% x3= 0 ;y3=-0;
% p1 = 0.083300; p2= 0.127889;
% vx1= p1; vy1= p2;
% vx2= p1;vy2= p2;
% vx3= -2*p1; vy3= -2*p2;


% BUTTERFLY I  需要高精度小步长
% http://three-body.ipb.ac.rs/sol.php?id=2
% x1= -1 ;y1= 0;
% x2= 1 ;y2= 0;
% x3= 0 ;y3=-0;
% p1 = 0.306893; p2= 0.125507;
% vx1= p1; vy1= p2;
% vx2= p1;vy2= p2;
% vx3= -2*p1; vy3= -2*p2;




planet1 = [x1,   vx1,   y1,  vy1];
planet2 = [x2,   vx2,   y2,  vy2];
planet3 = [x3,   vx3,   y3,  vy3];
t = [0,5.5];
ic = [planet1,planet2,planet3];
orbit(t,ic,t(2)*10000);

%%
function m1 = GetGlobal_m1()
% 用于获取第一个天体的质量
    m1 = 1;
end
function m2 = GetGlobal_m2()
% 用于获取第2个天体的质量
    m2 = 1;
end
function m3 = GetGlobal_m3()
% 用于获取第3个天体的质量
    m3 = 1;
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
    p = 100; % 每走100步就画图，这样节省空间，不用全部存储
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
        pause(0.001)
    end
    
    head1 = scatter(y(1,1),y(1,3),30,'red','filled');
    head2 = scatter(y(1,5),y(1,7),30,'blue','filled');
    head3 = scatter(y(1,9),y(1,11),30,'g','filled');
    legend([head1,head2,head3],'天体1','天体2','天体3',Location='best')
    axis auto
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