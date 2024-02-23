J = eye(4);
inv_J = eye(4);
u = zeros(4,1);
w = zeros(4,1);
count = 0;
while count < 100
    count = count + 1;
    % 以下利用Shermann-Morrison-Woodbury公式求逆
    inv_J_new = inv_J - inv_J*u*w'*inv_J./(1 + w'*inv_J*u);
    inv_J = inv_J_new;
    norm(inv_J*J-eye(4),'fro')
    
     u = 10^10*randn(4,1);
     w = 10^(-1)*randn(4,1);
    J2 = J + u*w'; % 更新J
    J = J2;
end