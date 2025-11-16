clc;
clear;
close all;

c1   = 1;
c2   = 1;
a12  = 0.005;
a21  = 0.0025;

% osobne równania
f1 = @(y1,y2)  c1*y1 - a12*y1*y2;         % dy1/dt
f2 = @(y1,y2) -c2*y2 + a21*y1*y2;         % dy2/dt

a = 0;
b = 25;

H = [0.005, 0.0025];      % wybrano mniejsze kroki, na potrzeby wykresu

y0 = [400; 80];

attr = {'Interpreter','latex'};

figure

for k = 1:2

    h = H(k);
    t = a:h:b;
    n = length(t);

    y = zeros(2, n);
    y(:,1) = y0;

    for i = 1:n-1
        y(1,i+1) = y(1,i) + h * f1(y(1,i), y(2,i));   % równanie 1
        y(2,i+1) = y(2,i) + h * f2(y(1,i), y(2,i));   % równanie 2
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
