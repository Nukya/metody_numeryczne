clc;
clear;
close all;

f = @(x) 1 ./(1 + 25*x.^2);

% funkcja interpolacji Lagrange'a
function s = lagrange(xp, yp, x)
    n = length(xp);
    s = 0;
    for k = 1:n
        p = 1;
        for i = 1:n
            if i ~= k
                p = p * (x - xp(i)) / (xp(k) - xp(i));
            end
        end
        s = s + yp(k) * p;
    end
end

ileWezlow = [3 5 6 8 9];
a = -1.5;
b = 1.5;

figure;
hold on;
grid on;

% funckja Czebyszewa
xd = linspace(a, b, 1000);
plot(xd, f(xd), 'k', 'LineWidth', 2);

legendEntries = {'$f(x) = \frac{1}{1 + 25x^2}$'};

blad = zeros(size(ileWezlow));

for i = 1:length(ileWezlow)
    n = ileWezlow(i);

    % węzły równoodległe
    xp = linspace(a, b, n);
    yp = f(xp);

    % wartości interpolacji
    yd = arrayfun(@(x) lagrange(xp, yp, x), xd);

    % średni błąd
    blad(i) = mean(abs(f(xd) - yd));

    % wykres interpolacji
    plot(xd, yd, 'LineWidth', 1.2);

    legendEntries{end+1} = sprintf('Interpolacja dla %d węzłów', n);
end

plot(xp, yp, 'ro', 'MarkerFaceColor', 'r');

title('Interpolacja metodą wielomianów Lagrangea dla równoodległych węzłów');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$f(x)$', 'Interpreter', 'latex');
legend(legendEntries, 'Interpreter', 'latex', 'Location', 'south');

bladRound = round(blad, 4, 'significant');

T = table(ileWezlow', bladRound', ...
    'VariableNames', {'Liczba węzłow', 'Średni bład interpolacji'});
disp(T);

% saveas(gcf, fullfile(pwd, 'zestaw_1_zadanie_4_2.png'));
