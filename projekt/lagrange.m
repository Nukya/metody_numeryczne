function s = lagrange(xp, yp, x)
%LAGRANGE Interpolacja wielomianowa metoda Lagrange'a.
%   s = LAGRANGE(xp, yp, x) zwraca wartosc wielomianu przechodzacego
%   przez punkty (xp, yp) obliczona w punkcie x. Funkcja pochodzi z
%   rozwiazan zestawu 1 (zadanie 4.2) i zostala wyodrebniona do osobnego
%   pliku w celu ponownego wykorzystania.

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
