function [t, y] = metodaEuleraUlepszona(fun, t0, tk, h, y0, params)
%METODAEULERAULEPSZONA   Wersja ulepszona (Heun/trapezy) z zestawów 1–3.
%   [t, y] = METODAEULERAULEPSZONA(fun, t0, tk, h, y0, params) rozwiązuje
%   układ y' = fun(t, y, params) w przedziale [t0, tk] krokiem h. Wektor y0
%   przechowuje warunki początkowe, a params jest strukturą parametrów
%   przekazywaną dalej do funkcji pochodnej.

    liczbaKrokow = floor((tk - t0) / h);
    czas = t0 + (0:liczbaKrokow)' * h;

    rozw = zeros(liczbaKrokow + 1, numel(y0));
    rozw(1, :) = y0(:).';

    for i = 1:liczbaKrokow
        k1 = fun(czas(i), rozw(i, :).', params);
        rozwProg = rozw(i, :) + h * k1.';
        k2 = fun(czas(i) + h, rozwProg.', params);
        rozw(i + 1, :) = rozw(i, :) + h / 2 * (k1.' + k2.');
    end

    if czas(end) < tk - eps
        czas(end + 1, 1) = tk;
        k1 = fun(czas(end - 1), rozw(end - 1, :).', params);
        rozwProg = rozw(end - 1, :) + (tk - czas(end - 1)) * k1.';
        k2 = fun(tk, rozwProg.', params);
        rozw(end, :) = rozw(end - 1, :) + (tk - czas(end - 1)) / 2 * (k1.' + k2.');
    end

    t = czas;
    y = rozw;
end
