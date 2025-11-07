clear;
close all;
clc;

% Dane

x_data = [-5 -4 -3 -2 -1 0 1 2 3 4 5 6];
y_data = [1.2 1.1 0.95 1.25 1.1 1.0 2 3.1 3.9 4.4 4.6 4.61];

a = -2;
b = 2;
liczbaProbek = 200;
x_fit = linspace(min(x_data), max(x_data), liczbaProbek);

degrees = 3:9;
results_polyfit = zeros(size(degrees));
results_mnk = zeros(size(degrees));

g_interp = @(x) interp1(x_data, y_data, x, 'linear', 'extrap');
f_interp = @(x) exp(-0.1 * x) .* (g_interp(x)).^2;
I_trapez = integral(f_interp, a, b);
disp(['Wartość całki (interpolacja liniowa): ', num2str(I_trapez)]);

% polyfit

figure;
plot(x_data, y_data, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Dane pomiarowe');
hold on;

colors = lines(length(degrees));

for i = 1:length(degrees)
    n = degrees(i);
    
    p = polyfit(x_data, y_data, n);
    
    g_polyfit = @(t) polyval(p, t);
    f = @(t) exp(-0.1*t) .* (g_polyfit(t)).^2;
    results_polyfit(i) = integral(f, a, b);
    
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Polyfit st. %d', n));
end

title('Aproksymacja Polyfit dla stopni 3–9');
xlabel('x');
ylabel('f(x)');
grid on;
legend('Location', 'northwest');
hold off;



% Metoda najmniejszych kwadratów
figure;
plot(x_data, y_data, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Dane pomiarowe');
hold on;

for i = 1:length(degrees)
    n = degrees(i);
    
    [wsp, ~] = aproksymacja_wielomianowa(x_data, y_data, n);
    
    g_mnk = @(t) polyval(wsp, t);
    f = @(t) exp(-0.1*t) .* (g_mnk(t)).^2;
    results_mnk(i) = integral(f, a, b);
    
    y_fit2 = polyval(wsp, x_fit);
    plot(x_fit, y_fit2, '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('MNK st. %d', n));
end

title('Aproksymacja MNK (własna metoda) dla stopni 3–9');
xlabel('x');
ylabel('f(x)');
grid on;
legend('Location', 'northwest');
hold off;

T = table(degrees', results_polyfit', results_mnk', ...
    'VariableNames', {'Stopień', 'Polyfit', 'MNK'});
disp(' ');
disp('===== Wartości całek dla obu metod =====');
disp(T);


function [wspolczynniki, wartosci_aproksymowane] = aproksymacja_wielomianowa(x, y, stopien)
    % Budowa macierzy Vandermonde'a (od najwyższego stopnia)
    y = y(:);
    A = zeros(length(x), stopien + 1);
    for i = 0:stopien
        A(:, stopien + 1 - i) = x.^i;
    end
    % Rozwiązanie równań normalnych
    wspolczynniki = (A' * A) \ (A' * y);
    wartosci_aproksymowane = A * wspolczynniki;
end
