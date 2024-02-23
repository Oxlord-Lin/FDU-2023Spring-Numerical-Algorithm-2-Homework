clear;
clc;
u = [1,1,1i];
v = [1,1i];
u = u(:);
v = v(:);
n1 = length(u)-1; % 最高次数
n2 = length(v)-1; % 最高次数
u = [u;zeros(n2,1)];
v = [v;zeros(n1,1)];
U = fft(u);
V = fft(v);
Y = U.*V;
y = ifft(Y)

