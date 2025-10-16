% Obliczanie wartości wielomianu interpolacyjnego Newtona

clc; 
clear; 
close all;

% Dane
x = [0.1 0.3 0.6 0.8];
y = [-1 1.2 1.0 1.5];
x_val = 0.55;

n = length(x);

F = zeros(n, n);     % stworzenie macierzy n x n
F(:,1) = y';         % uzupełnienie pierwszej kolumny macierzy wartosciami y

% Ilorazy n-1 rzędu
for j = 2:n
    for i = j:n
        F(i,j) = (F(i,j-1) - F(i-1,j-1)) / (x(i) - x(i-j+1));
    end
end

C = [x' F]

e = diag(F)';       % wyciągniecie współczynników z przekątnej

WN = e(1);
temp = 1;
for k = 2:n
    temp = temp * (x_val - x(k-1));
    WN = WN + e(k) * temp;
end

% przygotowanie danych do wykresu
xd = linspace(0, 1, 200);
yd = zeros(size(xd));

for j = 1:length(xd)
    temp2 = 1;
    yd(j) = e(1);
    for k = 2:n
        temp2 = temp2 * (xd(j) - x(k-1));
        yd(j) = yd(j) + e(k) * temp2;
    end
end

% wykres
figure;
plot(xd, yd, 'b-', 'LineWidth', 1.5); hold on;
plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot(x_val, WN, 'ks', 'MarkerFaceColor', 'g', 'MarkerSize', 8);

xlabel('x');
ylabel('WN(x)');
title('Interpolacja Newtona - Aproksymacja wielomianowa');
legend('Wielomian interpolacyjny', 'Punkty pomiarowe', sprintf('WN(%.3f)=%.3f', x_val, WN), 'Location', 'Best');
grid on;