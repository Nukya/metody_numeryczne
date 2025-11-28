function [t, y] = metodaEulera(fun, t0, tk, h, y0, params)
%METODAEULERA   Klasyczny schemat Eulera zgodny ze składnią z zestawów 1–3.
%   [t, y] = METODAEULERA(fun, t0, tk, h, y0, params) rozwiązuje układ
%   równań różniczkowych y' = fun(t, y, params) w przedziale [t0, tk]
%   krokiem h. Wektor y0 jest warunkiem początkowym, a params przenosi
%   dodatkowe parametry modelu (opcjonalne, może być pustą strukturą).

    liczbaKrokow = floor((tk - t0) / h);
    czas = t0 + (0:liczbaKrokow)' * h;

    rozw = zeros(liczbaKrokow + 1, numel(y0));
    rozw(1, :) = y0(:).';

    for i = 1:liczbaKrokow
        k1 = fun(czas(i), rozw(i, :).', params);
        rozw(i + 1, :) = rozw(i, :) + h * k1.';
    end

    if czas(end) < tk - eps
        czas(end + 1, 1) = tk;
        rozw(end, :) = rozw(end - 1, :) + (tk - czas(end - 1)) * fun(czas(end - 1), rozw(end - 1, :).', params).';
    end

    t = czas;
    y = rozw;
end
