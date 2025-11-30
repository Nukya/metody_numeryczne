function wart = wartosc_wielomianu(wsp, x)
% Oblicza wartość wielomianu w postaci Hornera dla współczynników wsp.
wart = 0;
for i = 1:length(wsp)
    wart = wart * x + wsp(i);
end
end
