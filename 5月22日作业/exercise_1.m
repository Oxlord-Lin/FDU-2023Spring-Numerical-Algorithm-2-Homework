ER = [];
H = [];
for n = 2:40
    h = pi/n;
    H = [H,h];
    s = 0;
    for k = 0:n-1
        s = s + sin( (k+1/3)*h ) + sin( (k+2/3)*h );
    end
    s = (h/2) *s;
    er = abs(2-s);
    ER = [ER,er];
end
plot(log(H),log(ER),'r-o',LineWidth=1);
xlabel('log(h)',FontSize=13,Interpreter='latex')
ylabel('log(Error)',FontSize=13,Interpreter='latex')
saveas(gcf,'exercise1-log(error)-against-log(h).png')