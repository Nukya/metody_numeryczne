function wspolczynniki = aproksymacja_wielomianowa(x, y, stopien)
%APROKSYMACJA_WIELOMIANOWA Dopasowanie wielomianu metoda MNK.
%   wspolczynniki = APROKSYMACJA_WIELOMIANOWA(x, y, stopien) zwraca
%   wspolczynniki wielomianu dopasowanego do danych (x, y) przy uzyciu
%   macierzy Vandermonde'a. Funkcja pochodzi z zestawu 1 (zad. 2) i
%   zostala przeniesiona do osobnego pliku do wielokrotnego uzycia.

    y = y(:);
    A = zeros(length(x), stopien + 1);
    for i = 0:stopien
        A(:, stopien + 1 - i) = x(:).^i;
    end
    wspolczynniki = (A' * A) \ (A' * y);
end
