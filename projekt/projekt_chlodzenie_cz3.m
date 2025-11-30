function projekt_chlodzenie_cz3
% Część 3 projektu – analiza sterowania procesem chłodzenia pręta przy
% zmiennym współczynniku h(DeltaT). Korzystamy wyłącznie z metod z
% "podręcznika" (interpolacja Lagrange'a, aproksymacja MNK, metoda Eulera
% ulepszona) oraz z funkcji przygotowanych w poprzednich częściach.

clc; close all;

%% Dane pomiarowe h(DeltaT) i dopasowanie MNK
[DeltaT_pom, h_pom] = pobierz_dane_pomiarowe();
[wsp_mnk, stopien_mnk] = konfiguruj_mnk(DeltaT_pom, h_pom);
fprintf('Wybrany stopień MNK: %d\n', stopien_mnk);
fprintf('Współczynniki MNK (od najwyższej potęgi):\n');
fprintf('  %.6f', wsp_mnk);
fprintf('\n\n');

h_model = @(Tb, Tw) h_mnk(Tb - Tw, wsp_mnk);

%% Parametry bazowe symulacji
param.A  = 0.0109;   % [m^2]
param.mb = 0.2;      % [kg]
param.mw = 2.5;      % [kg]
param.cb = 3.85;     % [kJ/(kg*K)]
param.cw = 4.1813;   % [kJ/(kg*K)]
param.Tb0 = 1200;    % [C]
param.Tw0 = 25;      % [C]
param.krok = 0.001;  % krok czasowy
param.tmax = 6;      % maks. czas symulacji

%% Zadanie 1: czas dojścia do zadanych progów temperatury
progi = [900 700 500 300 150];
[czasy_progi, przebiegi_ref] = czas_do_progow(progi, h_model, param);

fprintf('Czasy osiągnięcia progów temperatury (h = h_MNK):\n');
fprintf('%8s | %8s\n', 'Prog Tb', 't [s]');
fprintf(repmat('-',1,21)); fprintf('\n');
for i = 1:numel(progi)
    fprintf('%8.1f | %8.3f\n', progi(i), czasy_progi(i));
end

%% Zadanie 2: wpływ masy chłodziwa na dynamikę (planowanie sterowania)
mw_warianty = [1.5 2.5 3.5 5.0];
wyniki_masy = porownaj_masy(mw_warianty, h_model, param, 300);

fprintf('\nPorównanie czasu dojścia do 300C dla różnych mas chłodziwa:\n');
fprintf('%6s | %8s\n', 'mw [kg]', 't300 [s]');
fprintf(repmat('-',1,21)); fprintf('\n');
for i = 1:numel(mw_warianty)
    fprintf('%6.2f | %8.3f\n', mw_warianty(i), wyniki_masy.czasy(i));
end

%% Wizualizacje
rysuj_przebieg_progow(przebiegi_ref, progi);
rysuj_porownanie_mas(wyniki_masy, mw_warianty);

end

%% --- Zadanie 1 ---------------------------------------------------------
function [czasy, przebiegi] = czas_do_progow(progi, h_fun, param)
czasy = zeros(size(progi));
przebiegi = cell(size(progi));
for i = 1:numel(progi)
    prog = progi(i);
    [t, rozw] = symuluj_do_progu(prog, h_fun, param);
    czasy(i) = t(end);
    przebiegi{i} = struct('t', t, 'Tb', rozw(:,1), 'Tw', rozw(:,2), 'prog', prog);
end
end

function [t, rozw] = symuluj_do_progu(prog, h_fun, param)
odetmp = @(czas, x, par) f_zmienny_h(czas, x, par, h_fun);
[t, rozw] = metodaEuleraUlepszona(odetmp, 0, param.tmax, param.krok, [param.Tb0; param.Tw0], param);
ind = find(rozw(:,1) <= prog, 1, 'first');
if ~isempty(ind)
    t = t(1:ind);
    rozw = rozw(1:ind,:);
end
end

%% --- Zadanie 2 ---------------------------------------------------------
function wynik = porownaj_masy(mw_warianty, h_fun, param, prog)
liczba = numel(mw_warianty);
trajektorie = cell(liczba,1);
czasy = zeros(liczba,1);
for i = 1:liczba
    par = param;
    par.mw = mw_warianty(i);
    [t, rozw] = symuluj_do_progu(prog, h_fun, par);
    trajektorie{i} = struct('t', t, 'Tb', rozw(:,1), 'Tw', rozw(:,2), 'mw', par.mw);
    czasy(i) = t(end);
end
wynik.trajektorie = trajektorie;
wynik.czasy = czasy;
end

%% --- Modele h(DeltaT) --------------------------------------------------
function [wsp_mnk, stopien] = konfiguruj_mnk(DeltaT, h)
[~, wsp_mnk, stopien] = dopasuj_mnk(DeltaT, h);
end

function h_val = h_mnk(DeltaT, wsp_mnk)
DeltaT = DeltaT(:).';
h_val = zeros(size(DeltaT));
for i = 1:numel(DeltaT)
    h_val(i) = wartosc_wielomianu(wsp_mnk, DeltaT(i));
end
h_val = h_val(:);
end

%% --- Równania ODE ------------------------------------------------------
function dx = f_zmienny_h(~, x, param, h_fun)
Tb = x(1); Tw = x(2);
h = h_fun(Tb, Tw);
A = param.A; mb = param.mb; mw = param.mw; cb = param.cb; cw = param.cw;
dTb = -(h*A)/(mb*cb) * (Tb - Tw);
dTw =  +(h*A)/(mw*cw) * (Tb - Tw);
dx = [dTb; dTw];
end

%% --- Narzędzia MNK i wielomiany ---------------------------------------
function [err_struct, wsp, k_opt] = dopasuj_mnk(DeltaT, h)
stopnie = 2:8;
RMSE = zeros(size(stopnie));
wsp_store = cell(size(stopnie));
for i = 1:numel(stopnie)
    k = stopnie(i);
    wspolczynniki = aproksymacja_wielomianowa(DeltaT, h, k);
    aproks = h_mnk(DeltaT, wspolczynniki);
    RMSE(i) = sqrt(mean((aproks - h(:)).^2));
    wsp_store{i} = wspolczynniki;
end
improvement = (RMSE(1:end-1) - RMSE(2:end)) ./ RMSE(1:end-1);
idx = find(improvement < 0.01, 1, 'first');
if isempty(idx)
    [~, idx] = min(RMSE);
else
    idx = idx + 1;
end
k_opt = stopnie(idx);
wsp = wsp_store{idx};
err_struct = RMSE(idx);
end

function wart = wartosc_wielomianu(wsp, x)
wart = 0;
for i = 1:length(wsp)
    wart = wart * x + wsp(i);
end
end

%% --- Dane pomiarowe ----------------------------------------------------
function [DeltaT_pom, h_pom] = pobierz_dane_pomiarowe()
DeltaT_pom = [-1500 -1000 -300 -50 -1 1 20 50 200 400 1000 2000];
h_pom      = [178   176   168  161 160 160 160.2 161 165 168 174 179];
end

%% --- Wizualizacje ------------------------------------------------------
function rysuj_przebieg_progow(przebiegi, progi)
figure; hold on;
kolory = lines(numel(progi));
for i = 1:numel(progi)
    plot(przebiegi{i}.t, przebiegi{i}.Tb, 'LineWidth',1.3, 'Color', kolory(i,:), ...
        'DisplayName', sprintf('prog %.0fC', progi(i)));
end
xlabel('t [s]'); ylabel('T_b [C]');
title('Osiąganie zadanych progów temperatury');
legend('Location','best'); grid on;
end

function rysuj_porownanie_mas(wyniki, mw_warianty)
figure; hold on;
kolory = lines(numel(mw_warianty));
for i = 1:numel(mw_warianty)
    traj = wyniki.trajektorie{i};
    plot(traj.t, traj.Tb, 'LineWidth',1.3, 'Color', kolory(i,:), ...
        'DisplayName', sprintf('mw=%.1f kg', mw_warianty(i)));
end
xlabel('t [s]'); ylabel('T_b [C]');
title('Wpływ masy chłodziwa na czas dojścia do progu');
legend('Location','best'); grid on;
end

end
