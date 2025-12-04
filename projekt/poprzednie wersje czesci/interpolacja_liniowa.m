function yq = interpolacja_liniowa(x, y, xq)
% INTERPOLACJA_LINIOWA - własna implementacja interpolacji liniowej
%
%   x  - punkty znane (nierównoodległe), wektor 1×N
%   y  - wartości w punktach x, wektor 1×N
%   xq - punkty zapytania (wektor dowolnego rozmiaru)
%
%   yq - wynik interpolacji liniowej w punktach xq
%
% Metoda:
%   Dla każdego xq znajdowany jest przedział [x(i), x(i+1)],
%   a następnie stosowany jest wzór liniowy:
%
%   y = y_i + (y_{i+1}-y_i)/(x_{i+1}-x_i) * (xq - x_i)
%

    % upewniamy się, że x rośnie
    [x, idx] = sort(x);
    y = y(idx);

    yq = zeros(size(xq));

    for k = 1:numel(xq)
        xx = xq(k);

        % Jeśli xq jest poza zakresem — ekstrapolacja liniowa
        if xx <= x(1)
            yq(k) = y(1);
            continue;
        elseif xx >= x(end)
            yq(k) = y(end);
            continue;
        end

        % znajdź przedział
        for i = 1:length(x)-1
            if xx >= x(i) && xx <= x(i+1)
                % interpolacja liniowa
                t = (xx - x(i)) / (x(i+1) - x(i));
                yq(k) = y(i) * (1 - t) + y(i+1) * t;
                break;
            end
        end

    end

end
