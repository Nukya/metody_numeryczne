clear;
close all;
clc;

a = -1;
b = 2;
ilePrzedzialow = 6;
krok = (b-a)/ilePrzedzialow;

% Definicja funkcji podcałkowej
f  = @(x) sin(x) .* exp(x);

xi = linspace(a, b, ilePrzedzialow+1);
yi = f(xi);
srodekPrzedzialu = (xi(1:end-1) + xi(2:end))/2;
IProstokat = krok * sum(f(srodekPrzedzialu));

% Rysowanie funkcji podcałkowej wraz z dopasowanymi prostokątami

width = krok; % podstawa prostokąta
for i = 1:(length(xi)-1)
    posX = xi(i);
    height = f(srodekPrzedzialu(i));
    if height <=0
        posY = height;
        height = abs(height);

        hexStr = "#C492B1";
        color = hex2rgb(hexStr);
    else
        posY = 0;

        hexStr = "#BCE7FD";
        color = hex2rgb(hexStr);
    end
    rectangle('Position', [posX, posY, width, height], FaceColor = color);
end

hold on;

attr = {'Interpreter', 'latex'};
title('Metoda złożonych prostokątów');
subtitle('$I = \int_{-1}^{2} \sin(x)e^x \,dx$', 'Interpreter', 'latex');
plot(xi,yi,'b-','LineWidth', 2.0)
xlabel('$x$', attr{:});
ylabel('$y$', attr{:});

yline(0,'--r','LineWidth', 2.0)

hold off;
grid on; 

fprintf('Przybliżona wartość całki metodą złożonych prostokątów: %.4g\n', IProstokat);
