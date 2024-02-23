S = [];
H = [];
for n = 293:295
    h = 1/n;
    H = [H,h];
    s = 0;
    for k = 0:n-1
        s = s + exp((k+0.5) * h);
    end
    s = h * s;
    S = [S,s];
end
plot(H,S,'-o')