clear; 
clc; 
close all;

% Definicja funkcji p(x)
p = @(x) -0.1*x.^4 + 0.8*x.^3 - 0.6*x.^2 - 2*x + 1.5;

% Rozpatrywany przedział
a = -3;
b = 6;

% krok siatki
h = 0.1;

x_skan = a:h:b;
y_skan = p(x_skan);

pierwiastki = [];

ileIteracji = 30;

% Szukanie zmiany znaku
for i = 1:length(x_skan)-1
    if y_skan(i) * y_skan(i+1) < 0
        % Jeśli w przedziale zmienia się znak 
        % -> pierwiastek znajduje się w przedziale x_skan(i) and x_skan(i+1)
        left = x_skan(i);
        right = x_skan(i+1);

        % Bisekcja
        for k = 1:ileIteracji  % 30 iteracji
            mid = (left + right)/2;
            if p(left) * p(mid) <= 0
                right = mid;
            else
                left = mid;
            end
        end
        pierwiastki(end+1) = (left + right)/2;
    end
end

% Wyniki
disp('Znalezione pierwiastki funkcji p(x)');
disp(pierwiastki);

% Wykres
attr = {'Interpreter','latex'};

xd = linspace(a, b, 400);
yd = p(xd);
plot(xd, yd, 'b-', 'LineWidth', 1.5); hold on;
yline(0, '--k');
plot(pierwiastki, p(pierwiastki), 'ro', 'MarkerFaceColor','r');
xlabel('$x$', attr{:}); 
ylabel('$p(x)$', attr{:});
title('$p(x) = -0.1x.^4 + 0.8x.^3 - 0.6x.^2 - 2x + 1.5$', attr{:});
grid on;

for i = 1:length(pierwiastki)
    xr = pierwiastki(i);
    yr = p(xr);
    text(xr, yr-0.7, sprintf('(%.4g)', xr), ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold');
end