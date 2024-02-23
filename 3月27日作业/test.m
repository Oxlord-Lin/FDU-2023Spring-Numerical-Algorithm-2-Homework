Z = myfunc1(x,y,z,s,t)
function Z = myfunc1(x,y,z,s,t)
    point = [s,t]';
    n = length(x);
    denom = 0;
    for i = 1:n
        denom = denom + (1/norm([x(i),y(i)]'-point))
    end
    Z = 0;
    for i = 1:n
        Z = Z + z(i)*(norm([x(i),y(i)]'-point)/denom)
    end
end