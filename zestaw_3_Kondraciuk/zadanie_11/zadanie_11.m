clc;
clear;
close all;


% 1. Metoda złożona parabol Simpsona

m = 1000; % liczba przedzialow

%przedzial calkowania
a = -0.1;
b = 0.1;

h = (b - a) / m;
i = [0:1:m-2];

f = @(x) x .* (1 ./ sin((x + 0.5).^2)) .* exp(12 .* x); % funkcja podcałkowa


xd = linspace(a, b, m);
yd = f(xd);

suma = 0;
for iterator = 1:m-3
    suma = suma + h/3 * ( ...
                        yd(iterator) + ...
                        4 * yd(iterator+1) + ...
                        yd(iterator+2) ...
                        );
end

calka_simpson = suma;

% 2. Kwadratury Gaussa-Legendre'a n = 1,2,3,4

% Węzły Legendre'a na [-1,1]
t1 = [ -0.57735 ; 0.57735 ];
t2 = [ -0.774597 ; 0 ; 0.774597 ];
t3 = [ -0.861136 ; -0.339981 ; 0.339981 ; 0.861136 ];
t4 = [ -0.90618 ; -0.538469 ; 0 ; 0.538469 ; 0.90618 ];

% Wagi Legendre'a
A1 = [ 1 ; 1 ];
A2 = [ 5/9 ; 8/9 ; 5/9 ];
A3 = [ 0.347855 ; 0.652145 ; 0.652145 ; 0.347855 ];
A4 = [ 0.236927 ; 0.478629 ; 0.568889 ; 0.478629 ; 0.236927 ];

% funkcja anonimowa skalująca węzły na [a,b]
skaluj = @(t, a, b) (b - a)/2 .* t + (b + a)/2;

% funkcja anonimowa licząca całkę Gaussa dla danego zestawu wag i węzłów
gauss_calc = @(Ai, xi, a, b) (b - a)/2 * sum(Ai .* f(xi));

% Skalowanie węzłów
x1 = skaluj(t1, a, b);
x2 = skaluj(t2, a, b);
x3 = skaluj(t3, a, b);
x4 = skaluj(t4, a, b);

% Obliczenia Gaussa
calka_1 = gauss_calc(A1, x1, a, b);
calka_2 = gauss_calc(A2, x2, a, b);
calka_3 = gauss_calc(A3, x3, a, b);
calka_4 = gauss_calc(A4, x4, a, b);

% 3. Wizualizacja funkcji i aproksymacji

x_plot = linspace(a, b, 2000);
y_plot = f(x_plot);

% Wielomiany interpolacyjne (1–4 stopnia)
p1 = polyfit(x1, f(x1), 1); y1 = polyval(p1, x_plot);
p2 = polyfit(x2, f(x2), 2); y2 = polyval(p2, x_plot);
p3 = polyfit(x3, f(x3), 3); y3 = polyval(p3, x_plot);
p4 = polyfit(x4, f(x4), 4); y4 = polyval(p4, x_plot);

figure;
plot(x_plot, y_plot, 'k', 'LineWidth', 1.4); 
hold on;
plot(x_plot, y1, '--', 'LineWidth', 1.2);
plot(x_plot, y2, '--', 'LineWidth', 1.2);
plot(x_plot, y3, '--', 'LineWidth', 1.2);
plot(x_plot, y4, '--', 'LineWidth', 1.2);

attr = {'Interpreter','latex'};

legend(["$f(x)$", "wiel. st. 1", "wiel. st. 2", "wiel. st. 3", "wiel. st. 4"], ...
       'Location','best', attr{:});

title('Porownanie funkcji $f(x)$ i aproksymacji wielomianowych $n=1,2,3,4$ stopnia', ...
      attr{:});
xlabel("$x$", attr{:});
ylabel("$f(x)$", attr{:});
grid on;

% 4. Tabela wyników

Metoda = ["Węzły Legendre'a"; "Węzły Legendre'a"; "Węzły Legendre'a"; "Węzły Legendre'a"; ...
          "Metoda złożona parabol"; "Wynik Wolfram"];
Rzad = [1;2;3;4;1000;0];
Wartosc = [calka_1; calka_2; calka_3; calka_4; calka_simpson; 0.0237525];

T = table(Metoda, Rzad, Wartosc);
disp(T);

% saveas(gca, "zadanie_11.png");
