clear
close all
clc

f = @(x) 1 ./ (1 + exp(-x));

xi = -5:1:5;
yi = f(xi);

plot(xi, yi)
grid on
xlabel('x')
ylabel('f(x)')
title('f(x) = 1 / (1 + e^{-x})')

%  Budowa macierzy (bez x^2) 
X = [ones(size(xi(:))), xi(:), xi(:).^3, xi(:).^4];

%  Obliczenie współczynników wielomianu w sensie najmniejszych kwadratów
a = X \ yi(:);
a0 = a(1); a1 = a(2); a3 = a(3); a4 = a(4);


p = @(x) a0+a1.*x + a3.*x.^3 + a4.*x.^4;

xd = linspace(-5, 5, 1000);
yd = f(xd);
pd = p(xd);

% =================

diff = yd - pd;               % wektor różnic
MAE = mean(abs(diff));        % średni błąd bezwzględny 
%                               -> zsumuj wektor różnic i podziel prze liczbę punktów
RMSE = sqrt(mean(diff.^2));   % pierwiastek z błędu średniokwadratowego

%Wyniki 
fprintf('Współczynniki aproksymacji:\n');
fprintf('a0 = %.6f\n', a0);
fprintf('a1 = %.6f\n', a1);
fprintf('a3 = %.6f\n', a3);
fprintf('a4 = %.6f\n', a4);
fprintf('\nŚredni błąd bezwzględny = %.6f\n', MAE);
fprintf('Pierwiastek z błędu średniokwadratowego (RMSE) = %.6f\n', RMSE);

%Wykres 
attr = {'Interpreter', 'latex'};
figure('Color','w');
plot(xd, yd, 'LineWidth', 1.5, 'DisplayName', '$f(x) = \frac{1}{(1+e^{-x})}$');
hold on;
plot(xd, pd, '--', 'LineWidth', 1.5, 'DisplayName', '$p(x) = a_0 + a_1x + a_3x^3 + a_4x^4$');
plot(xi, yi, 'ro', 'DisplayName', 'węzły aproksymacji');
grid on;
xlabel('$x$', attr{:});
ylabel('$y$', attr{:});
title('Aproksymacja wielomianowa 3-go stopnia');
legend('Location','best', attr{:});
set(gca, 'FontSize', 12);

% saveas(gcf,'zadanie_7.png')

