clc;
clear;
close all;

f = @(x) 1/5 + 1/20*x + 1./(1+exp(-2*(x-1))) - 1./(1+exp(-3*x));

x = linspace(-5,5,100);
y = f(x);

A = [];
B = [];
R = [];

for k = 1:length(x)-1
    if y(k)*y(k+1) <= 0
        A = [A x(k)];
        B = [B x(k+1)];
    end
end

for m = 1:length(A)
    a = A(m);
    b = B(m);
    fa = f(a);
    fb = f(b);
    xm = (a*fb - b*fa)/(fb - fa);
    fm = f(xm);
    i = 1;
    while i < 1000 && abs(fm) >= 1e-6
        if fa*fm <= 0
            b = xm;
            fb = fm;
        else
            a = xm;
            fa = fm;
        end
        xm = (a*fb - b*fa)/(fb - fa);
        fm = f(xm);
        i = i+1;
    end
    R = [R xm];
end

attr = {'Interpreter','latex'};

plot(x, y)
hold on
yline(0)
plot(R, f(R), 'o')
xlabel('$x$', attr{:})
ylabel('$f(x)$', attr{:})
hold off

saveas(gca, "zadanie_14.png");