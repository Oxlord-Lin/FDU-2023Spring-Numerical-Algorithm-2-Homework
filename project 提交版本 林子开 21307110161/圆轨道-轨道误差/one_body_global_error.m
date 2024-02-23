ns = [1000,2000,3000,5000,8000,1e4,2e4,3e4,5e4,8e4,1e5];
hs = [];
ex_trap_error = [];
SV_error = [];
RK4_error = [];
for k = 1:length(ns)
    n = ns(k);
    [z,h] = orbit([0,10],[0,2*pi,1,0],n,2);
    hs = [hs,h];
    ex_trap_error = [ex_trap_error,norm([0,1]-z)];
    [z,~] = orbit([0,10],[0,2*pi,1,0],n,3);
    SV_error = [SV_error,norm([0,1]-z)];
    [z,~] = orbit([0,10],[0,2*pi,1,0],n,4);
    RK4_error = [RK4_error,norm([0,1]-z)];
end
p1 = loglog(hs,ex_trap_error,'-o',LineWidth=1); hold on
p2 = loglog(hs,SV_error,'--h',LineWidth=1);
p3 = loglog(hs,RK4_error,'-.^',LineWidth=1);
legend([p1,p2,p3],'显式梯形法','Störmer-Verlet法','RK4',Location='best')
xlabel('log(h)',FontSize=13)
ylabel('log(global error)',FontSize=13)
title('显式梯形法，Störmer-Verlet法，RK4的全局误差随步长h的变化（双对数图）')
saveas(gcf,'三种方法的全局误差随步长h的变化（双对数图）.png')

[z,h] = orbit([0,1],[0,2*pi,1,0],10000,3);

% 2：显式梯形（好），3：Störmer-Verlet（好），4：4阶龙格-库塔（好）

function g = GetGlobal_g()
% 获得万有引力常数
    g = 1; % 归一化了，这使得在数量上有m2*g=m2，便于计算
end

function m2 = GetGlobal_m2()
% 获得中心天体质量
    m2 = 4*pi^2; % 在这个部分这样设置，使得周期为1
end

function [z,h]=orbit(inter,ic,n,method)
% input如下
% 时间区间 inter
% 初始条件 ic=[x0,vx0,y0,vy0],分别表示初始x坐标，初始x方向速度，初始y坐标，初始y方向速度
% 步数n
% 调用单步方法,用method指定
% 2：显式梯形（好），3：Störmer-Verlet（好），4：4阶龙格-库塔（好）
% output:z是最后所在坐标，h是步长
    h = (inter(2)-inter(1))/n; % 步长
    x0 = ic(1); vx0 = ic(2); y0 = ic(3); vy0 = ic(4); % 初始条件
    y = [x0,vx0,y0,vy0];
    t(1) = inter(1);
%     % 显式欧拉法
%     if method == 1
%         for k = 1:n
%             t(k+1) = t(k) + h;
%             y(k+1,:) = eulerstep(t(k),y(k,:),h);
%         end

    % 显式梯形法
    if method == 2
        for k = 1:n
            t_new = t + h;
            y = trapstep(t,y,h);
            t = t_new;
        end
        z = y(end,:);
        z = [z(1),z(3)];

    % Störmer-Verlet法
    elseif method == 3
        v = [vx0,vy0];
        r = [x0,y0];
        v = SV_v_step(t,r(1,:),v(1,:),0.5*h); % 第一次，速度走"半步"
        for k = 1:n
            t_new = t + h;
            r = r + h*v;
            v = SV_v_step(t,r,v,h);
            t = t_new;
        end
        z = r(end,:);

    % 4阶龙格-库塔法
    elseif method == 4
        for k = 1:n
            t_new = t + h;
            y = RK4step(t,y,h);
            t = t_new;
        end
        z = y(end,:);
        z = [z(1),z(3)];
%     % 隐式欧拉法
%     elseif method == 5
%         for k = 1:n
%             t(k+1) = t(k) + h;
%             y_init = eulerstep(t(k),y(k,:),h); % 先用显式欧拉走一步作为牛顿法初值
%             y(k+1,:) = implicit_Euler_Newton(t(k),y(k,:),y_init,h); % 隐式欧拉法，需要用牛顿法求根
%         end
%     % 隐式梯形法
%     elseif method == 6
%         Y = y(1,:);
%         x0 = Y(1); vx0 = Y(2); y0 = Y(3); vy0 = Y(4);
%         m2 = 3; g = 1; 
%         hmg = h*m2*g; dist3=(x0^2+y0^2)^(3/2);
%         f = [h*vx0, hmg*(-x0)/dist3, h*vy0, hmg*(-y0)/dist3]'; % 初始时左端点处的dY/dt
%         for k = 1:n
%             t(k+1) = t(k) + h;
%             y_init = trapstep(t(k),y(k,:),h); % 先用显式梯形走一步作为牛顿法初值
%             [f_new,y(k+1,:)] = implicit_trap_Newton(t(k),f,y(k,:),y_init,h); % 隐式梯形法，需要用牛顿法求根,f=dY/dt
%             f = f_new;
%         end
    end
end

% function y = eulerstep(t,x,h)
% % 显示欧拉法
% % input:x是当前状态，t是当前时间，h是步长
% % output：y是下一刻的状态
%     y = x + h*ydot(t,x);  % ydot是x处的导数
% end

function y = trapstep(t,x,h)
% 显式梯形法
% input:x是当前状态，t是当前时间，h是步长
% output：y是下一刻的状态
    z1 = ydot(t,x);
    xtemp = x+h*z1;
    z2 = ydot(t+h,xtemp);
    y = x + h*(z1 + z2)/2;
end

function v_next = SV_v_step(t,r,v,h)
% Störmer-Verlet法，本函数用于更新v
% input:r是k步位置，v是(k-1/2)步的速度，t是当前时间，h是步长
% output：v_next是(k+1/2)步的速度
    m2 = GetGlobal_m2();
    g = GetGlobal_g();
    mg2 = m2*g; px2 = 0; py2 = 0;
    px1 = r(1); py1 = r(2); vx1 = v(1); vy1 = v(2);
    dist = sqrt((px2-px1)^2 + (py2-py1)^2);
    v_next = zeros(1,2);
    v_next(1) = vx1 + h * (mg2*(px2-px1)) / dist^3;
    v_next(2) = vy1 + h * (mg2*(py2-py1)) / dist^3;
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

function z = ydot(t,x)
% input:x是当前状态，t是当前时间，h是步长
% output：ydot是x处的导数，即dy/dt
    m2 = GetGlobal_m2();
    g = GetGlobal_g();
    mg2 = m2*g; px2 = 0; py2 = 0;
    px1 = x(1); py1 = x(3); vx1 = x(2); vy1 = x(4);
    dist = sqrt((px2-px1)^2 + (py2-py1)^2);
    z = zeros(1,4);
    z(1) = vx1;
    z(2) = (mg2*(px2-px1))/dist^3;
    z(3) = vy1;
    z(4) = (mg2*(py2-py1))/dist^3;
end

% function Y2 = implicit_Euler_Newton(t,Y0,Y,h)
% % input:Y0是上一步的状态，Y是牛顿法的初值，t是当前时间，h是步长
% % output：Y2是牛顿法求出的根
%     m2 = 3; g = 1; 
%     hmg = h*m2*g;
%     Y = Y(:);
%     Y0 = Y0(:);
%     Fnorm = 1;
%     t
%     while Fnorm > 1e-14
%         x = Y(1); vx = Y(2); y = Y(3); vy = Y(4); dist5=(x^2+y^2)^(5/2);dist3=(x^2+y^2)^(3/2);
%         F = Y - Y0 - [h*vx, hmg*(-x)/dist3, h*vy, hmg*(-y)/dist3]';
%         Fnorm = norm(F);
%         J = [1,                    -h,   0,                    0;
%             hmg*(y^2-2*x^2)/dist5,  1,   -3*hmg*x*y/dist5,    0;
%             0,                     0,   1,                     -h;
%             -3*hmg*x*y/dist5,      0,   hmg*(x^2-2*y^2)/dist5,  1];
%         Y = Y - (J\F);
%     end
%     Y2 = Y;
% end
% 
% function [f_new,Y2] = implicit_trap_Newton(t,f,Y0,Y,h)
% % input:Y0是左端点，Y是牛顿法的初值，t是当前时间，f是左端点的导数，h是步长
% % output：Y2是牛顿法求出的根，f_new=dY2/dt,是右端点的导数值
%     m2 = 3; g = 1; 
%     hmg = h*m2*g;
%     Y = Y(:);
%     Y0 = Y0(:);
%     Fnorm = 1;
%     t
%     while Fnorm > 1e-14
%         x = Y(1); vx = Y(2); y = Y(3); vy = Y(4); dist5=(x^2+y^2)^(5/2);dist3=(x^2+y^2)^(3/2);
%         F = Y - Y0 - 0.5*h*f -0.5*[h*vx, hmg*(-x)/dist3, h*vy, hmg*(-y)/dist3]';
%         Fnorm = norm(F);
%         J = [1,                        -h/2,    0,                          0;
%             0.5*hmg*(y^2-2*x^2)/dist5,  1,      -0.5*3*hmg*x*y/dist5,       0;
%             0,                          0,      1,                         -h/2;
%             -0.5*3*hmg*x*y/dist5,       0,      0.5*hmg*(x^2-2*y^2)/dist5,  1];
%         Y = Y - (J\F);
%     end
%     Y2 = Y;
%     x = Y(1); vx = Y(2); y = Y(3); vy = Y(4) ;dist3=(x^2+y^2)^(3/2);
%     f_new = [h*vx, hmg*(-x)/dist3, h*vy, hmg*(-y)/dist3]';
% end