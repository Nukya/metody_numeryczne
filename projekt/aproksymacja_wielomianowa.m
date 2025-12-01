function [wspolczynniki, wartosci_aproksymowane] = aproksymacja_wielomianowa(x, y, stopien)
    y = y(:);
    A = zeros(length(x), stopien + 1);

    % Macierz Vandermonde'a (od najwyższego stopnia)
    for i = 0:stopien
        A(:, stopien + 1 - i) = x.^i;
    end

    % Rozwiązanie równań normalnych
    wspolczynniki = (A' * A) \ (A' * y);

    wartosci_aproksymowane = A * wspolczynniki;
end
