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
