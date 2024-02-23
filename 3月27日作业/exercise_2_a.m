data = [ 1, 1 
    0, 2 
    -1, 1 
    0, 0 
    1, -1 
    0, -2 
    -1, -1 ];

x = data(:,1);
y = data(:,2);
[n,~] = size(data);
t = zeros(n,1);
for i = 2:n
    t(i) = norm(data(i,:)-data(i-1,:))+t(i-1);
end

tt = 0:.0001:t(end);
xt = spline(t,x,tt);
yt = spline(t,y,tt);
plot(xt,yt,LineWidth=1)
hold on
plot(data(:,1),data(:,2),'*')
