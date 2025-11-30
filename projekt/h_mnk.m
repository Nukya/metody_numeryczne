function h_val = h_mnk(DeltaT, wsp_mnk)
% Wylicza wartości h(DeltaT) dla wielomianu dopasowanego metodą MNK.
DeltaT = DeltaT(:).';
h_val = zeros(size(DeltaT));
for i = 1:numel(DeltaT)
    h_val(i) = wartosc_wielomianu(wsp_mnk, DeltaT(i));
end
h_val = h_val(:);
end
