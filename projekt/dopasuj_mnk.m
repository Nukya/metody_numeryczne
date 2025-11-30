function [err_struct, wsp, k_opt] = dopasuj_mnk(DeltaT, h)
% Dopasowuje wielomian metodą najmniejszych kwadratów i wybiera stopień.
stopnie = 2:8;
RMSE = zeros(size(stopnie));
wsp_store = cell(size(stopnie));
for i = 1:numel(stopnie)
    k = stopnie(i);
    wspolczynniki = aproksymacja_wielomianowa(DeltaT, h, k); % zestaw 1 zad. 2
    aproks = h_mnk(DeltaT, wspolczynniki);
    RMSE(i) = sqrt(mean((aproks - h(:)).^2));
    wsp_store{i} = wspolczynniki;
end
improvement = (RMSE(1:end-1) - RMSE(2:end)) ./ RMSE(1:end-1);
idx = find(improvement < 0.01, 1, 'first');
if isempty(idx)
    [~, idx] = min(RMSE);
else
    idx = idx + 1; % stopień odpowiada kolejnemu elementowi
end
k_opt = stopnie(idx);
wsp = wsp_store{idx};
err_struct = RMSE(idx);
end
