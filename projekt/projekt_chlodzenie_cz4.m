function projekt_chlodzenie_cz4
% Czwarta część projektu – proste sterowanie czasem chłodzenia oraz
% utrzymanie zadanej temperatury przy zmiennym współczynniku h(DeltaT).
% Wykorzystujemy wyłącznie narzędzia z "podręcznika" (interpolacja
% Lagrange'a, aproksymacja MNK, metoda Eulera ulepszona) oraz funkcje
% przygotowane w poprzednich częściach projektu.

clc; close all;

%% Dane pomiarowe h(DeltaT) i dopasowanie MNK
[DeltaT_pom, h_pom] = pobierz_dane_pomiarowe();
[coef_mnk, stopien_mnk] = konfiguruj_mnk(DeltaT_pom, h_pom);
fprintf('Wybrany stopień MNK: %d\n', stopien_mnk);
fprintf('Współczynniki MNK (od najwyższej potęgi):\n');
fprintf('  %.6f', coef_mnk);
fprintf('\n\n');

h_model = @(Tb, Tw) h_mnk(Tb - Tw, coef_mnk);

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

%% Zadanie 1: dobór masy chłodziwa do wymaganego czasu dojścia do progu
prog_Tb = 300;
t_docel = 3.5; % [s]
[mw_opt, traj_opt] = znajdz_mase_do_czasu(prog_Tb, t_docel, h_model, param);

fprintf('=== Dobór masy dla czasu dojścia ===\n');
fprintf('Docelowy czas osiągnięcia %.0fC: %.3f s\n', prog_Tb, t_docel);
fprintf('Znaleziono masę chłodziwa: %.3f kg\n', mw_opt);
fprintf('Rzeczywisty czas dojścia: %.3f s\n\n', traj_opt.t(end));

rysuj_trajektorie_masy(traj_opt, prog_Tb);

%% Zadanie 2: sterowanie dwustanowe utrzymujące temperaturę blisko zadanego poziomu
T_zad = 250;    % [C]
pasmowosc = 5;  % histereza +/- 5C
t_sym = 5;      % czas symulacji
wynik_ctrl = sterowanie_dwustanowe(T_zad, pasmowosc, h_model, param, t_sym);

fprintf('=== Sterowanie dwustanowe ===\n');
fprintf('Temperatura zadana: %.1fC, pasmo: +/- %.1fC\n', T_zad, pasmowosc);
fprintf('Czas zejścia poniżej %.1fC: %.3f s\n', T_zad, wynik_ctrl.t_osiagniecia);
fprintf('Maksymalne przeregulowanie: %.3fC\n', wynik_ctrl.przeregulowanie);
fprintf('Średni błąd po osiągnięciu zadania: %.3fC\n', wynik_ctrl.sredni_blad);

rysuj_sterowanie_dwustanowe(wynik_ctrl, T_zad, pasmowosc);

end

%% --- Zadanie 1: odwrotne zagadnienie masy ------------------------------
function [mw_opt, traj] = znajdz_mase_do_czasu(prog, t_docel, h_fun, param)
% proste wyszukiwanie dwumianowe w zakresie 1–8 kg
lewy = 1.0; prawy = 8.0;
traj = symuluj_do_progu(prog, h_fun, aktualizuj_mw(param, (lewy+prawy)/2));
for iter = 1:12
    srodek = (lewy + prawy) / 2;
    par = aktualizuj_mw(param, srodek);
    traj_local = symuluj_do_progu(prog, h_fun, par);
    czas = traj_local.t(end);
    if czas > t_docel
        % za wolno – zwiększamy masę chłodziwa (większe chłodzenie)
        lewy = srodek;
    else
        prawy = srodek;
        traj = traj_local;
    end
end
mw_opt = (lewy + prawy) / 2;
traj.param = aktualizuj_mw(param, mw_opt);
end

function traj = symuluj_do_progu(prog, h_fun, param)
ode = @(czas, x, par) f_zmienny_h(czas, x, par, h_fun);
[t, rozw] = metodaEuleraUlepszona(ode, 0, param.tmax, param.krok, [param.Tb0; param.Tw0], param);
ind = find(rozw(:,1) <= prog, 1, 'first');
if isempty(ind)
    ind = length(t);
end
traj.t = t(1:ind);
traj.Tb = rozw(1:ind,1);
traj.Tw = rozw(1:ind,2);
end

function param = aktualizuj_mw(param, mw)
param.mw = mw;
end

%% --- Zadanie 2: sterowanie dwustanowe ----------------------------------
function wynik = sterowanie_dwustanowe(T_zad, pasmo, h_fun, param, tmax)
Tb = param.Tb0; Tw = param.Tw0;
krok = param.krok;
N = floor(tmax / krok);

% Przygotowanie tablic wynikowych
t = zeros(N+1,1); Tb_hist = zeros(N+1,1); Tw_hist = zeros(N+1,1); u_hist = zeros(N+1,1);
Tb_hist(1) = Tb; Tw_hist(1) = Tw; u_hist(1) = 1; % startujemy z chłodzeniem wł.

for k = 1:N
    % Decyzja sterownika z histerezą
    if Tb > T_zad + pasmo
        u_hist(k) = 1;
    elseif Tb < T_zad - pasmo
        u_hist(k) = 0;
    else
        u_hist(k) = u_hist(k-1+ (k==1)); % utrzymaj poprzedni stan (dla k=1 zostaje 1)
    end

    h_eff = h_z_uwzglednieniem_sterowania(h_fun, Tb, Tw, u_hist(k));
    k1 = rhs(Tb, Tw, h_eff, param);
    Tb_tmp = Tb + krok * k1(1);
    Tw_tmp = Tw + krok * k1(2);
    h_eff2 = h_z_uwzglednieniem_sterowania(h_fun, Tb_tmp, Tw_tmp, u_hist(k));
    k2 = rhs(Tb_tmp, Tw_tmp, h_eff2, param);

    Tb = Tb + (krok/2) * (k1(1) + k2(1));
    Tw = Tw + (krok/2) * (k1(2) + k2(2));

    t(k+1) = t(k) + krok;
    Tb_hist(k+1) = Tb;
    Tw_hist(k+1) = Tw;
end

% ostatnia decyzja dla czytelności wykresu
u_hist(end) = u_hist(end-1);

% Analiza wyników
idx_osiagniecia = find(Tb_hist <= T_zad, 1, 'first');
if isempty(idx_osiagniecia)
    idx_osiagniecia = length(t);
end
blad_po = Tb_hist(idx_osiagniecia:end) - T_zad;

wynik.t = t;
wynik.Tb = Tb_hist;
wynik.Tw = Tw_hist;
wynik.u = u_hist;
wynik.t_osiagniecia = t(idx_osiagniecia);
wynik.przeregulowanie = max(Tb_hist) - T_zad;
wynik.sredni_blad = mean(abs(blad_po));
end

function h_val = h_z_uwzglednieniem_sterowania(h_fun, Tb, Tw, u)
% u=1: pełne chłodzenie, u=0: podtrzymanie (10% strumienia)
h_on = h_fun(Tb, Tw);
h_off = 0.1 * h_on + 10; % niewielkie chłodzenie podtrzymujące
h_val = u * h_on + (1 - u) * h_off;
end

function dx = rhs(Tb, Tw, h, param)
A = param.A; mb = param.mb; mw = param.mw; cb = param.cb; cw = param.cw;
dTb = -(h*A)/(mb*cb) * (Tb - Tw);
dTw =  +(h*A)/(mw*cw) * (Tb - Tw);
dx = [dTb; dTw];
end

%% --- Modele h(DeltaT) ---------------------------------------------------
function [wsp_mnk, stopien] = konfiguruj_mnk(DeltaT, h)
[~, wsp_mnk, stopien] = dopasuj_mnk(DeltaT, h);
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

%% --- Wizualizacje ------------------------------------------------------
function rysuj_trajektorie_masy(traj, prog)
figure; hold on;
plot(traj.t, traj.Tb, 'LineWidth',1.3,'DisplayName','T_b');
plot(traj.t, traj.Tw, 'LineWidth',1.3,'DisplayName','T_w');
plot(traj.t, prog*ones(size(traj.t)), 'k--','DisplayName',sprintf('prog %.0fC', prog));
xlabel('t [s]'); ylabel('T [C]');
title(sprintf('Dobór masy chłodziwa: mw = %.3f kg', traj.param.mw));
legend('Location','best'); grid on;
end

function rysuj_sterowanie_dwustanowe(wynik, T_zad, pasmo)
figure;
subplot(2,1,1); hold on;
plot(wynik.t, wynik.Tb, 'LineWidth',1.3,'DisplayName','T_b');
plot(wynik.t, wynik.Tw, 'LineWidth',1.3,'DisplayName','T_w');
plot(wynik.t, (T_zad+pasmo)*ones(size(wynik.t)), 'k--','DisplayName','pasmo');
plot(wynik.t, (T_zad-pasmo)*ones(size(wynik.t)), 'k--','HandleVisibility','off');
xlabel('t [s]'); ylabel('T [C]'); grid on; legend('Location','best');
title('Sterowanie dwustanowe – przebiegi temperatur');

subplot(2,1,2);
stairs(wynik.t, wynik.u, 'LineWidth',1.3);
ylabel('u'); xlabel('t [s]'); ylim([-0.1 1.1]); grid on;
title('Sygnał sterujący (1 – chłodzenie włączone)');
end

end
