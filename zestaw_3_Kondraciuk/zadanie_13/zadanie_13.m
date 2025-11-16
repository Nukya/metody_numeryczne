clc;
clear;
close all;

c1   = 1;
c2   = 1;
a12  = 0.005;
a21  = 0.0025;

f1 = @(y1,y2) c1*y1 - a12*y1*y2;
f2 = @(y1,y2) -c2*y2 + a21*y1*y2;

a = 0;
b = 25;

H = [0.005, 0.0025];
y0 = [400; 80];

attr = {'Interpreter','latex'};

figure

for k = 1:2

    h = H(k);
    t = a:h:b;
    n = length(t);

    y1(1) = y0(1);
    y2(1) = y0(2);

    for i = 1:n-1
        k1_1 = f1(y1(i), y2(i));
        k1_2 = f2(y1(i), y2(i));

        y1_mid = y1(i) + h/2 * k1_1;
        y2_mid = y2(i) + h/2 * k1_2;

        k2_1 = f1(y1_mid, y2_mid);
        k2_2 = f2(y1_mid, y2_mid);

        y1(i+1) = y1(i) + h * k2_1;
        y2(i+1) = y2(i) + h * k2_2;
    end

    subplot(2,1,k)
    hold on
    plot(t, y1, 'b', 'DisplayName', '$y_1(t)$')
    plot(t, y2, 'r', 'DisplayName', '$y_2(t)$')

    title(['$h = ', num2str(h), '$'], attr{:})
    xlabel('$t$', attr{:})
    ylabel('$y(t)$', attr{:})
    legend(attr{:})
    grid on

end

saveas(gca, "zadanie_13.png");
