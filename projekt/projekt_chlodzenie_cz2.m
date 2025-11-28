function projekt_chlodzenie_cz2
% Rozszerzona symulacja chlodzenia preta w oleju z nieliniowym wspolczynnikiem h(DeltaT)
% W pliku zebrano trzy modele h(DeltaT): interpolacja wielomianowa, aproksymacja MNK i splajn kubiczny.
% Kazdy model zostaje wykorzystany w rozwiazaniu ODE opisujacym zmiane temperatur pręta i oleju.

clc; close all;

%% Dane pomiarowe h(DeltaT)
[DeltaT_pom, h_pom] = pobierz_dane_pomiarowe();

%% Parametry modelu fizycznego i symulacji
param.A  = 0.0109;   % [m^2]
param.mb = 0.2;      % [kg]
param.mw = 2.5;      % [kg]
param.cb = 3.85;     % [kJ/(kg*K)]
param.cw = 4.1813;   % [kJ/(kg*K)]
param.h_const = 160; % referencyjny wspolczynnik h
param.Tb0 = 1200;    % [C]
param.Tw0 = 25;      % [C]
param.tspan = [0 5];

%% Przygotowanie danych do modeli
[DeltaT_rowne, h_rowne] = generuj_rowne_odstepy(DeltaT_pom, h_pom, 200);
konfiguruj_mnk(DeltaT_pom, h_pom);        % zapis koeficjentow MNK w pamieci
konfiguruj_spline(DeltaT_rowne, h_rowne); % zapis danych do splajnu w pamieci

%% Porownanie h(DeltaT)
DeltaT_siatka = linspace(min(DeltaT_pom), max(DeltaT_pom), 400);
h_poly    = h_interp_poly(DeltaT_siatka);
h_mnk_fit = h_mnk(DeltaT_siatka);
h_spl     = h_spline(DeltaT_siatka);

%% Symulacje ODE dla kazdego modelu h
[t_const, rozw_const] = ode45(@(t,x) f_const(t,x,param), param.tspan, [param.Tb0; param.Tw0]);
[t_poly,  rozw_poly]  = ode45(@(t,x) f_interp_poly(t,x,param), param.tspan, [param.Tb0; param.Tw0]);
[t_mnk,   rozw_mnk]   = ode45(@(t,x) f_mnk(t,x,param), param.tspan, [param.Tb0; param.Tw0]);
[t_spl,   rozw_spl]   = ode45(@(t,x) f_spline(t,x,param), param.tspan, [param.Tb0; param.Tw0]);

%% Obliczanie bledow h(DeltaT) wzgledem danych pomiarowych
err_h_poly = oblicz_bledy(h_pom, h_interp_poly(DeltaT_pom));
err_h_mnk  = oblicz_bledy(h_pom, h_mnk(DeltaT_pom));
err_h_spl  = oblicz_bledy(h_pom, h_spline(DeltaT_pom));

%% Obliczanie bledow Tb(t) wzgledem modelu referencyjnego h = const
[t_ref, Tb_ref] = deal(t_const, rozw_const(:,1));
Tb_poly = interp1(t_poly, rozw_poly(:,1), t_ref, 'pchip');
Tb_mnk  = interp1(t_mnk,  rozw_mnk(:,1),  t_ref, 'pchip');
Tb_spl  = interp1(t_spl,  rozw_spl(:,1),  t_ref, 'pchip');
err_Tb_poly = oblicz_bledy(Tb_ref, Tb_poly);
err_Tb_mnk  = oblicz_bledy(Tb_ref, Tb_mnk);
err_Tb_spl  = oblicz_bledy(Tb_ref, Tb_spl);

%% Wizualizacje
rysuj_h(DeltaT_pom, h_pom, DeltaT_siatka, h_poly, h_mnk_fit, h_spl);
rysuj_temperatury(t_const, rozw_const, t_poly, rozw_poly, t_mnk, rozw_mnk, t_spl, rozw_spl, param.h_const);
rysuj_bledy_Tb(t_ref, Tb_poly, Tb_mnk, Tb_spl, Tb_ref);

%% Raport tekstowy
fprintf('\nTabela bledow h(DeltaT) wzgledem pomiarow:\n');
drukuj_tabele_bledow({'Interpolacja poly','MNK','Splajn'}, [err_h_poly; err_h_mnk; err_h_spl]);

fprintf('\nTabela bledow Tb(t) wzgledem modelu h=const:\n');
drukuj_tabele_bledow({'Interpolacja poly','MNK','Splajn'}, [err_Tb_poly; err_Tb_mnk; err_Tb_spl]);

fprintf('\nWnioski:\n');
interpretuj_wyniki(err_h_poly, err_h_mnk, err_h_spl, err_Tb_poly, err_Tb_mnk, err_Tb_spl);

end

%% --- Modele h(DeltaT) -------------------------------------------------
function h_val = h_interp_poly(DeltaT)
persistent DeltaT_pom h_pom;
if isempty(DeltaT_pom)
    [DeltaT_pom, h_pom] = pobierz_dane_pomiarowe();
end
% interpolacja wielomianowa metoda Lagrange'a (z zestaw 1 zad. 4.2)
h_val = arrayfun(@(dt) lagrange(DeltaT_pom, h_pom, dt), DeltaT);
end

function konfiguruj_mnk(DeltaT, h)
persistent wsp_mnk stopien;
[~, wspolczynniki, k] = dopasuj_mnk(DeltaT, h);
wsp_mnk = wspolczynniki;
stopien = k;
fprintf('Wybrany stopien MNK: %d\n', stopien);
fprintf('Wspolczynniki MNK (od najwyzszej potegi):\n');
fprintf('  %.6f', wsp_mnk);
fprintf('\n');
assignin('caller','mnk_stopien',stopien);
assignin('caller','mnk_wspolczynniki',wsp_mnk);
end

function h_val = h_mnk(DeltaT)
persistent wsp_mnk stopien;
if isempty(wsp_mnk)
    [DeltaT_pom, h_pom] = pobierz_dane_pomiarowe();
    [~, wsp_mnk, stopien] = dopasuj_mnk(DeltaT_pom, h_pom);
end
h_val = polyval(wsp_mnk, DeltaT);
end

function konfiguruj_spline(DeltaT_rowne, h_rowne)
persistent pp_spline;
pp_spline = spline(DeltaT_rowne, h_rowne);
assignin('caller','pp_spline',pp_spline);
end

function h_val = h_spline(DeltaT)
persistent pp_spline;
if isempty(pp_spline)
    [DeltaT_pom, h_pom] = pobierz_dane_pomiarowe();
    [DeltaT_rowne, h_rowne] = generuj_rowne_odstepy(DeltaT_pom, h_pom, 200);
    pp_spline = spline(DeltaT_rowne, h_rowne);
end
h_val = ppval(pp_spline, DeltaT);
end

%% --- Funkcje ODE ------------------------------------------------------
function dx = f_const(~, x, param)
Tb = x(1); Tw = x(2);
h = param.h_const;
dx = ode_rhs(Tb, Tw, h, param);
end

function dx = f_interp_poly(~, x, param)
Tb = x(1); Tw = x(2);
h = h_interp_poly(Tb - Tw);
dx = ode_rhs(Tb, Tw, h, param);
end

function dx = f_mnk(~, x, param)
Tb = x(1); Tw = x(2);
h = h_mnk(Tb - Tw);
dx = ode_rhs(Tb, Tw, h, param);
end

function dx = f_spline(~, x, param)
Tb = x(1); Tw = x(2);
h = h_spline(Tb - Tw);
dx = ode_rhs(Tb, Tw, h, param);
end

function dx = ode_rhs(Tb, Tw, h, param)
A = param.A;
mb = param.mb; mw = param.mw;
cb = param.cb; cw = param.cw;
dTb = -(h*A)/(mb*cb) * (Tb - Tw);
dTw =  +(h*A)/(mw*cw) * (Tb - Tw);
dx = [dTb; dTw];
end

%% --- Narzedzia matematyczne ------------------------------------------
function [err_struct, wsp, k_opt] = dopasuj_mnk(DeltaT, h)
stopnie = 2:8;
RMSE = zeros(size(stopnie));
wsp_store = cell(size(stopnie));
for i = 1:numel(stopnie)
    k = stopnie(i);
    wspolczynniki = aproksymacja_wielomianowa(DeltaT, h, k); % zestaw 1 zad. 2
    aproks = polyval(wspolczynniki, DeltaT);
    RMSE(i) = sqrt(mean((aproks - h).^2));
    wsp_store{i} = wspolczynniki;
end
improvement = (RMSE(1:end-1) - RMSE(2:end)) ./ RMSE(1:end-1);
idx = find(improvement < 0.01, 1, 'first');
if isempty(idx)
    [~, idx] = min(RMSE);
else
    idx = idx + 1; % stopien odpowiada kolejnemu elementowi
end
k_opt = stopnie(idx);
wsp = wsp_store{idx};
err_struct = RMSE(idx);
end

function [DeltaT_rowne, h_rowne] = generuj_rowne_odstepy(DeltaT_pom, h_pom, N)
DeltaT_rowne = linspace(min(DeltaT_pom), max(DeltaT_pom), N);
h_rowne = interp1(DeltaT_pom, h_pom, DeltaT_rowne, 'linear');
end

function bledy = oblicz_bledy(ref, approx)
roznica = approx - ref;
bledy.MSE = mean(roznica.^2);
bledy.max = max(abs(roznica));
bledy.rel_mean = mean(abs(roznica)) / mean(abs(ref));
end

%% --- Wizualizacje i raporty ------------------------------------------
function rysuj_h(DeltaT_pom, h_pom, DeltaT_siatka, h_poly, h_mnk_fit, h_spl)
figure; hold on;
scatter(DeltaT_pom, h_pom, 40, 'k', 'filled', 'DisplayName','pomiary');
plot(DeltaT_siatka, h_poly, 'LineWidth',1.4,'DisplayName','interpolacja poly');
plot(DeltaT_siatka, h_mnk_fit,'LineWidth',1.4,'DisplayName','MNK');
plot(DeltaT_siatka, h_spl, 'LineWidth',1.4,'DisplayName','splajn');
plot(DeltaT_siatka, 160*ones(size(DeltaT_siatka)), 'k--','DisplayName','h const=160');
xlabel('\DeltaT [C]'); ylabel('h [W/(m^2*K)]');
title('Modele h(\DeltaT)'); grid on; legend('Location','best');
end

function rysuj_temperatury(t_const, rozw_const, t_poly, rozw_poly, t_mnk, rozw_mnk, t_spl, rozw_spl, h_const)
figure;
subplot(2,1,1); hold on;
plot(t_const, rozw_const(:,1),'k--','LineWidth',1.2,'DisplayName','h const');
plot(t_poly,  rozw_poly(:,1), 'LineWidth',1.2,'DisplayName','interpolacja poly');
plot(t_mnk,   rozw_mnk(:,1),  'LineWidth',1.2,'DisplayName','MNK');
plot(t_spl,   rozw_spl(:,1),  'LineWidth',1.2,'DisplayName','splajn');
title(sprintf('T_b(t) dla roznych h(\DeltaT) (h const = %g)', h_const));
xlabel('t [s]'); ylabel('T_b [C]'); grid on; legend('Location','best');

subplot(2,1,2); hold on;
plot(t_const, rozw_const(:,2),'k--','LineWidth',1.2,'DisplayName','h const');
plot(t_poly,  rozw_poly(:,2), 'LineWidth',1.2,'DisplayName','interpolacja poly');
plot(t_mnk,   rozw_mnk(:,2),  'LineWidth',1.2,'DisplayName','MNK');
plot(t_spl,   rozw_spl(:,2),  'LineWidth',1.2,'DisplayName','splajn');
title('T_w(t) dla roznych modeli h');
xlabel('t [s]'); ylabel('T_w [C]'); grid on; legend('Location','best');
end

function rysuj_bledy_Tb(t_ref, Tb_poly, Tb_mnk, Tb_spl, Tb_ref)
figure; hold on;
plot(t_ref, Tb_poly - Tb_ref, 'LineWidth',1.2,'DisplayName','interpolacja poly');
plot(t_ref, Tb_mnk  - Tb_ref, 'LineWidth',1.2,'DisplayName','MNK');
plot(t_ref, Tb_spl  - Tb_ref, 'LineWidth',1.2,'DisplayName','splajn');
plot(t_ref, zeros(size(t_ref)), 'k--','DisplayName','odniesienie');
title('Blad T_b(t) wzgledem modelu h=const');
xlabel('t [s]'); ylabel('\Delta T_b [C]'); grid on; legend('Location','best');
end

function drukuj_tabele_bledow(nazwy, bledy_struct)
M = numel(nazwy);
fprintf('%-20s | %10s | %10s | %12s\n', 'Model', 'MSE', 'max', 'rel_mean');
fprintf(repmat('-',1,62)); fprintf('\n');
for i = 1:M
    b = bledy_struct(i);
    fprintf('%-20s | %10.4f | %10.4f | %12.4f\n', nazwy{i}, b.MSE, b.max, b.rel_mean);
end
end

function interpretuj_wyniki(err_h_poly, err_h_mnk, err_h_spl, err_Tb_poly, err_Tb_mnk, err_Tb_spl)
err_h = [err_h_poly.MSE, err_h_mnk.MSE, err_h_spl.MSE];
err_T = [err_Tb_poly.MSE, err_Tb_mnk.MSE, err_Tb_spl.MSE];
[~, idx_h] = min(err_h);
[~, idx_T] = min(err_T);
model_names = {'interpolacja wielomianowa','MNK','splajn kubiczny'};
fprintf('Najlepsze dopasowanie h(DeltaT): %s\n', model_names{idx_h});
fprintf('Najmniejszy blad przebiegu T_b(t): %s\n', model_names{idx_T});
if idx_h == idx_T
    fprintf('Spójny wybór modelu: %s jest rekomendowany do symulacji.\n', model_names{idx_h});
else
    fprintf('Roznice w kryteriach wskazuja na potrzebe dodatkowej analizy (h: %s, T_b: %s).\n', model_names{idx_h}, model_names{idx_T});
end
end

%% --- Dane pomocnicze --------------------------------------------------
function [DeltaT_pom, h_pom] = pobierz_dane_pomiarowe()
DeltaT_pom = [-1500 -1000 -300 -50 -1 1 20 50 200 400 1000 2000];
h_pom      = [178   176   168  161 160 160 160.2 161 165 168 174 179];
end

end
