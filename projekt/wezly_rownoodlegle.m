function hq = wezly_rownoodlegle(dT, h, dTq)
% WEZLY_ROWNOODLEGLE - interpolacja Lagrange'a 1 stopnia (odcinkowa)
%
%  dT  - punkty znane (rosnące)
%  h   - wartości w punktach dT
%  dTq - punkty, w których liczymy h(dTq)
%
%  hq  - wartość interpolacji L1(dTq)

    % upewniamy się, że węzły są rosnące (opcjonalne)
    [dT, idx] = sort(dT);
    h = h(idx);

    m = numel(dTq);
    hq = zeros(size(dTq));

    for k = 1:m
        xq = dTq(k);

        % sprawdzamy zakres danych
        if xq < dT(1) || xq > dT(end)
            error('Punkt %.4f poza zakresem danych.', xq);
        end

        % szukamy jednego jedynego przedziału
        for i = 1:length(dT)-1

            if xq >= dT(i) && xq <= dT(i+1)

                % dwa punkty tworzące przedział
                x0 = dT(i);
                x1 = dT(i+1);
                y0 = h(i);
                y1 = h(i+1);

                % Lagrange 1 stopnia = interpolacja liniowa
                t = (xq - x0) / (x1 - x0);
                hq(k) = y0*(1 - t) + y1*t;

                break;  % mamy wynik, wychodzimy
            end

        end   
    end
end
