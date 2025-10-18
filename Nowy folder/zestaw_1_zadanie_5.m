clc;
clear;
close all;

xi = [0.2 0.4 0.6];
yi = [-1.79 -2.26 -1.59];

A = [4 2 0;
    1 4 1;
    0 2 4];

a = 0.2;
b = 0.6;
n = length(xi)-1;   % liczba węzłów
h = (b-a)/n;        %skok
alpha = 4.02;       % pochodna w punkcie 0.2
beta = -2.62;       % pochodna w punkcie 0.6

d = [yi(1) + alpha * h / 3 
    yi(2) 
    yi(3) - beta * h / 3];

c = A \ d;

c0 = c(1);
c1 = c(2);
c2 = c(3);

c_m1 = c1 - h / 3 * alpha;
c3 = c1 + h / 3 * beta;


% xd = linspace(0.2, 0.6, 100)
% phi = phi_fun(x5wezlow, h, xd); 
% S3 = @(x) c_m1 .* phi(-1, x) + c0 .* phi(0, x) + c1 .* phi(1, x) + c2 .* phi(2, x) + c3 .* phi(3, x);


 a = 0.2;
 b = 0.6;
 xd = linspace(a,b,100);
 xk = xi
 
 yd = zeros(length(xd),1);
 
 % używamy hold on, poniewaz chcemy dorysowac kolejne przebiegi
 % na tym samym rysunku
 hold on
 for i = -1:5
 %for i = 0:4
 %for i = 1:3
     % obliczamy wspolrzedna wezla srodkowego i-tej funkcji
     % bazowej
     xi = xk(1) + (i)*h;
     % narysujmy pionowa linie aby zaznaczyc srodek funkcji bazowej
     xline(xi, '--', {sprintf('x_{%d}',i)})  

      for j = 1:length(xd)
     yd(j) = phi_fun(xi, h, xd(j));
 end
plot(xd,yd)

hold off   

% plotToPrint = plot(xi,yi,"o");
% hold on
% t = 0:0.01:0.2;
% plot (t,S3(t))
% t = 0.2:0.01:0.4;
% plot (t,S3(t))
% t = 0.4:0.01:0.6;
% plot (t,S3(t))
% t = 0.6:0.01:0.8;
% plot (t,S3(t))
% grid on;
% 
% hold off
% 
% title 'xyz z5'
% xlabel 'x'
% ylabel 'y'

% saveas(plotToPrint,[pwd 'zestaw_1_zadanie_5.png'])

function y = phi_fun(xi, h, x)
    % xi - punkt środkowy funkcji bazowej
    % h - odległość między węzłami
    % x - wspolrzedna x, dla ktorej obliczamy wartosc

    if (x < xi - 2*h) || (x > xi + 2*h)
        y = 0;
    elseif x < xi - h
        y = (x - (xi - 2*h))^3;
    elseif x < xi
        y = (x - (xi - 2*h))^3 - 4*(x - (xi - h))^3;
    elseif x < xi + h
        y = ((xi + 2*h) - x)^3 - 4*((xi + h) - x)^3;
    else
        y = ((xi + 2*h) - x)^3;
    end
    y = y / h^3;
end