clc
clear
close all;
clc;
clear;

f = @(y1,y2) [ 1*y1 - 0.005*y1*y2;
              -1*y2 + 0.0025*y1*y2 ];

a = 0;
b = 25;

H = [0.01, 0.005];
y0 = [400; 80];

attr = {'Interpreter','latex'};

figure

for k = 1:2

    h = H(k);
    n = (b-a)/h;
    t = a:h:b;

    y = zeros(2, n+1);
    y(:,1) = y0;

    for i = 1:n
        y(:,i+1) = y(:,i) + h * f(y(1,i), y(2,i));
    end

    subplot(2,1,k)
    hold on
    plot(t, y(1,:), 'b', 'DisplayName', '$y_1(t)$')
    plot(t, y(2,:), 'r', 'DisplayName', '$y_2(t)$')

    title(['$h = ', num2str(h), '$'], attr{:})
    xlabel('$t$', attr{:})
    ylabel('$y(t)$', attr{:})
    legend(attr{:})
    grid on

end

% saveas(gca, "zadanie_12.png");
