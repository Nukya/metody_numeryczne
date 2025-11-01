clear;
clc;
close all;

% Definicja funkcji p(x)
p = @(x) -0.1*x.^4 + 0.8*x.^3 - 0.6*x.^2 - 2*x + 1.5;

% Pochodna funkcji p'(x)
dp = @(x) -0.4*x.^3 + 2.4*x.^2 - 1.2*x - 2;

% Rozpatrywany przedział
a = -3;
b = 6;

% Parametry iteracyjne
tol = 1e-8;      % dokładność
ileIteracji = 30;    % maksymalna liczba iteracji
h = 0.1;         % krok siatki do wykrywania zmian znaku

x_skan = a:h:b;
y_skan = p(x_skan);

pierwiastki = [];

% Szukanie miejsc zmiany znaku, aby wybrać dobre punkty startowe
for i = 1:length(x_skan)-1
    if y_skan(i) * y_skan(i+1) < 0
        % początkowy punkt startowy x0
        x0 = x_skan(i);
        for k = 1:ileIteracji
            if abs(dp(x0)) < 1e-12
                
                % Unikamy dzielenia przez zero
                break;
            end
            x1 = x0 - p(x0)/dp(x0);
            if abs(x1 - x0) < tol
                break;
            end
            x0 = x1;
        end
        pierwiastki(end+1) = x1;
    end
end

% Wyniki
disp('Znalezione pierwiastki funkcji p(x):');
disp(pierwiastki');

% Wykres funkcji i pierwiastków
attr = {'Interpreter','latex'};
xd = linspace(a, b, 400);
yd = p(xd);

figure;
plot(xd, yd, 'b-', 'LineWidth', 1.5); hold on;
yline(0, '--k');
plot(pierwiastki, p(pierwiastki), 'ro', 'MarkerFaceColor', 'r');
xlabel('$x$', attr{:});
ylabel('$p(x)$', attr{:});
title('$p(x) = -0.1x^4 + 0.8x^3 - 0.6x^2 - 2x + 1.5$', attr{:});
subtitle('$p''(x) = -0.4x^3 + 2.4x^2 - 1.2x - 2$', attr{:});
grid on;

for i = 1:length(pierwiastki)
    xr = pierwiastki(i);
    yr = p(xr);
    text(xr, yr - 0.7, sprintf('(%.4g)', xr), ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold');
end
