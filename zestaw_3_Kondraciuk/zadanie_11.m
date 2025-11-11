clc;
clear;
close all;

%===========================
% metoda złożona parabol Simpsona

m = 1000; % liczba przedzialow

%przedzial calkowania
a = -0.1;
b = 0.1;

h = (b - a) / m;
i = [0:1:m-2];

f = @(x) x .* 1 / (sin((x + 0.5).^2))*exp(12 .* x); % funkcja podcałkowa

xd = linspace(-0.1, 0.1, m);
yd = f(xd);

suma = 0;
for iterator = 1:m-3
    suma = suma + h/3 * ( ...
                        yd(iterator) + ...
                        4 * yd(iterator+1) + ...
                        yd(iterator+2) ...
                        );
end

disp(suma);