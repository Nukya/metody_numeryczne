close all;
clc;
clear;

xi = [0.2 0.4 0.6];
yi = [-1.79 -2.26 -1.59];

A = [4 2 0;
    1 4 1;
    0 2 4];

a = xi(1);
b = xi(end);
n = length(xi)-1;
h = (b-a)/n;
alpha = 4.02; % pochodna w punkcie 0.2
beta = -2.62; % pochodna w punkcie 0.6

d = [yi(1) + h / 3 * alpha 
    yi(2) 
    yi(3) - h / 3 * beta];

wspC = A \ d;   % wspolczynniki c od 0 do 2

c0 = wspC(1);
c1 = wspC(2);
c2 = wspC(3);
c_m1 = c1 - h/3*alpha;
c3 = c1 + h/3*beta;

wspolczynnikiC = [c_m1 c0 c1 c2 c3]; % wspolczynniki c od -1 do 3 
xiRozszerzone = [xi(1) - h, xi(1), xi(2), xi(3), xi(end) + h];

phiLista = zeros(size(xiRozszerzone));

for i = 1:length(xiRozszerzone)
    phiLista(i) = phi(xiRozszerzone(i), h, 0.23);
end

S3 = sum(wspolczynnikiC .* phiLista); % wynik dla S3(0.23)
fprintf('S3(x) = %.4f\n', S3);


% funkcja sklejona
xd = linspace(a - h, b + h, 300);
S = zeros(size(xd));
for k = 1:length(xiRozszerzone)
    S = S + wspolczynnikiC(k) * arrayfun(@(x) phi(xiRozszerzone(k), h, x), xd);
end


figure;
hold on; 
grid on;

% funkcje bazowe
for k = 1:length(xiRozszerzone)
    phi_all(k, :) = arrayfun(@(x) phi(xiRozszerzone(k), h, x), xd);
    plot(xd, phi_all(k,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('$\\phi_{%d}(x)$', k-2));
end

xlabel('$x$', 'Interpreter','latex');
ylabel('$\phi_i(x)$', 'Interpreter','latex');
title('Funkcje bazowe $\phi_i(x)$', 'Interpreter','latex');
legend('Interpreter','latex', 'Location','best');


plot(xd, S, 'b', 'LineWidth', 2, 'DisplayName', '$S_3(x)$');
plot(xi, yi, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'węzły');
plot(0.23, S3, 'ro', 'MarkerFaceColor', 'g', 'DisplayName', '$x=0.23$');

xlabel('$x$', 'Interpreter','latex');
ylabel('$S_3(x)$', 'Interpreter','latex');
title('Funkcja interpolujaca $S_3(x)$', 'Interpreter','latex');
legend show;
legend('Interpreter','latex','Location','northeastoutside');
grid on;






function wZakresie = czyWZakresie(x, a, b)
% czyWZakresie - zwraca true, jeśli x mieści się w zakresie [a,b]
    wZakresie = (x >= a) & (x <= b);
end


function result = phi(xi, h, x)
% xi - środek funkcji bazowej (czyli punkt x_i)
% h - odległość między węzłami
% x - punkt, w którym liczymy wartość

    % Wyznaczenie sąsiednich punktów węzłowych
    xi_m2 = xi - 2*h;
    xi_m1 = xi - h;
    xi_p1 = xi + h;
    xi_p2 = xi + 2*h;

    if czyWZakresie(x, xi_m2, xi_m1)
        y = (x - xi_m2)^3;
    elseif czyWZakresie(x, xi_m1, xi)
        y = (x - xi_m2)^3 - 4*(x - xi_m1)^3;
    elseif czyWZakresie(x, xi, xi_p1)
        y = (xi_p2 - x)^3 - 4*(xi_p1 - x)^3;
    elseif czyWZakresie(x, xi_p1, xi_p2)
        y = (xi_p2 - x)^3;
    else
        y = 0;
    end

    %Φ_i(x) = y / h^3
    result = y / h^3;
end



saveas(gcf, fullfile(pwd, 'zestaw_1_zadanie_5.png'));
