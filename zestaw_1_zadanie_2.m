clear;
close all;
clc;

% Dane
x_data = [-5 -4 -3 -2 -1 0 1 2 3 4 5 6];
y_data = [1.2 1.1 0.95 1.25 1.1 1.0 2 3.1 3.9 4.4 4.6 4.61];

% obliczenie wartości całki złożoną metodą trapezów

g = @(x) interp1(x_data, y_data, x, 'linear', 'extrap');

a = -2;
b = 2;
liczbaProbek = 100;
h = (b - a) / liczbaProbek;
x = linspace(a, b, liczbaProbek+1);
y = exp(-0.1 * x) .* (g(x)).^2;

I_trapez = h * (0.5*y(1) + sum(y(2:end-1)) + 0.5*y(end));
disp(['Wartość całki (złożona metoda trapezów): ', num2str(I_trapez)]);

% obliczenie wartości całki za pomocą polyfit

degrees = 3:9;

results = zeros(size(degrees));

for i = 1:length(degrees)
    n = degrees(i);
    
    % Aproksymacja wielomianowa
    p = polyfit(x_data, y_data, n);
    
    % Funkcja aproksymująca
    g = @(t) polyval(p, t);
    
    % Definicja funkcji podcałkowej
    f = @(t) exp(-0.1*t) .* (g(t)).^2;
    
    % Obliczenie całki numerycznie w zakresie (-2, 2)
    results(i) = integral(f, -2, 2);
end

x_fit = linspace(min(x_data), max(x_data), liczbaProbek); 
y_fit = polyval(p, x_fit);

figure;
plot(x_fit, y_fit, 'b-', 'LineWidth', 2); 
hold on;

% obliczenie wartości całki za pomocą własnej funkcji
% (aproksymacja wielomianowa metodą najmniejszych kwadratów)

for i = 1:length(degrees)
    n = degrees(i);

[wsp, y_apr] = aproksymacja_wielomianowa(x_data, y_data, n);

% Funkcja aproksymująca
g = @(t) polyval(wsp, t);

% Definicja funkcji podcałkowej
f = @(t) exp(-0.1*t) .* (g(t)).^2;

% Obliczenie całki numerycznie w zakresie (-2, 2)
wyniki(i) = integral(f, -2, 2);
end

y_fit2 = polyval(wsp, x_fit);


plot(x_data, y_data, 'o');
plot(x_fit, y_fit2);

xlim([-6 6.5])

legend('Aproksymacja', 'Dane');
xlabel('x', 'Interpreter', 'latex');
ylabel('f(x)', 'Interpreter', 'latex');
title('Aproksymacja');
hold off;

T = table(degrees', results', wyniki',...
    'VariableNames', {'Stopień', 'Polyfit', 'MNK'});
disp(T);


function [wspolczynniki, wartosci_aproksymowane] = aproksymacja_wielomianowa(x, y, stopien)
% APROKSYMACJA_WIELOMIANOWA Aproksymacja danych (x, y) wielomianem danego stopnia
% przy użyciu metody najmniejszych kwadratów.
%
% Wejście:
%   x - wektor argumentów (np. pomiarów)
%   y - wektor wartości funkcji (np. wyników pomiarów)
%   stopien - stopień wielomianu aproksymującego
%
% Wyjście:
%   wspolczynniki - współczynniki wielomianu (od najwyższego stopnia)
%   wartosci_aproksymowane - wartości wielomianu w punktach x

    y = y(:);

    % Budowa macierzy Vandermonde'a
    A = zeros(length(x), stopien + 1);
    for i = 0:stopien
        A(:, stopien + 1 - i) = x.^i;
    end

    % Rozwiązanie układu równań normalnych
    wspolczynniki = (A' * A) \ (A' * y);

    % Obliczenie wartości aproksymowanych
    wartosci_aproksymowane = A * wspolczynniki;
end



% saveas(gcf, 'D:\repos\metody_numeryczne\zestaw_1_zadanie_2.png')