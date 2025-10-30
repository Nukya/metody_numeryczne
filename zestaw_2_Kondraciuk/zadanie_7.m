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

%  Budowa macierzy projektującej (bez x^2) 
X = [ones(size(xi(:))), xi(:), xi(:).^3, xi(:).^4];

%  Obliczenie współczynników wielomianu w sensie najmniejszych kwadratów
a = X \ yi(:);
a0 = a(1); a1 = a(2); a3 = a(3); a4 = a(4);


p = @(x) a0+a1.*x + a3.*x.^3 + a4.*x.^4;

xd = linspace(-5, 5, 1000);
yd = f(xd);
pd = p(xd);

% ---------------------------------------------------

diff = yd - pd;               % wektor różnic
MAE = mean(abs(diff));        % średni błąd bezwzględny 
%                               -> zsumuj wektor różnic i podziel prze liczbę punktów
RMSE = sqrt(mean(diff.^2));   % pierwiastek z błędu średniokwadratowego

% --- Wyniki ---
fprintf('Współczynniki aproksymacji:\n');
fprintf('a0 = %.6f\n', a0);
fprintf('a1 = %.6f\n', a1);
fprintf('a3 = %.6f\n', a3);
fprintf('a4 = %.6f\n', a4);
fprintf('\nŚredni błąd bezwzględny = %.6f\n', MAE);
fprintf('Pierwiastek z błędu średniokwadratowego (RMSE) = %.6f\n', RMSE);

% --- Wykres ---
figure('Color','w');
plot(xd, yd, 'LineWidth', 1.5, 'DisplayName', 'f(x) = 1/(1+e^{-x})');
hold on;
plot(xd, pd, '--', 'LineWidth', 1.5, 'DisplayName', 'p(x) = a0 + a1x + a3x^3 + a4x^4');
plot(xi, yi, 'ro', 'DisplayName', 'węzły aproksymacji');
grid on;
xlabel('x'); ylabel('y');
title('Aproksymacja funkcji sigmoidalnej wielomianem bez składnika x^2');
legend('Location','best');
set(gca, 'FontSize', 12);