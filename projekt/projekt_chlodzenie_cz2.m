function projekt_chlodzenie_cz2
% Rozszerzona symulacja chłodzenia pręta z ograniczeniem do metod z
% podrecznika (zestawy 1–3). W kodzie stosujemy jedynie własne procedury
% interpolacji wielomianowej, aproksymacji MNK i całkowania metodą Eulera
% ulepszonego.

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
param.krok = 0.001;  % krok czasowy dla wszystkich modeli

%% Przygotowanie modeli h(DeltaT)
[wsp_mnk, stopien_mnk] = konfiguruj_mnk(DeltaT_pom, h_pom);
fprintf('Wybrany stopień MNK: %d\n', stopien_mnk);
fprintf('Współczynniki MNK (od najwyższej potęgi):\n');
fprintf('  %.6f', wsp_mnk);
fprintf('\n');

DeltaT_siatka = linspace(min(DeltaT_pom), max(DeltaT_pom), 400);
h_poly = h_interp_poly(DeltaT_siatka, DeltaT_pom, h_pom);
h_mnk_fit = h_mnk(DeltaT_siatka, wsp_mnk);
h_lin = h_interp_liniowa(DeltaT_siatka, DeltaT_pom, h_pom);

%% Symulacje ODE dla każdego modelu h (metoda Eulera ulepszona)
h_const_fun = @(Tb, Tw) param.h_const;
h_poly_fun  = @(Tb, Tw) h_interp_poly(Tb - Tw, DeltaT_pom, h_pom);
h_mnk_fun   = @(Tb, Tw) h_mnk(Tb - Tw, wsp_mnk);
h_lin_fun   = @(Tb, Tw) h_interp_liniowa(Tb - Tw, DeltaT_pom, h_pom);

[t_const, rozw_const] = symuluj_model(h_const_fun, param);
[t_poly,  rozw_poly]  = symuluj_model(h_poly_fun,  param);
[t_mnk,   rozw_mnk]   = symuluj_model(h_mnk_fun,   param);
[t_lin,   rozw_lin]   = symuluj_model(h_lin_fun,   param);

%% Obliczanie błędów h(DeltaT) względem danych pomiarowych
err_h_poly = oblicz_bledy(h_pom, h_interp_poly(DeltaT_pom, DeltaT_pom, h_pom));
err_h_mnk  = oblicz_bledy(h_pom, h_mnk(DeltaT_pom, wsp_mnk));
err_h_lin  = oblicz_bledy(h_pom, h_interp_liniowa(DeltaT_pom, DeltaT_pom, h_pom));

%% Obliczanie błędów Tb(t) względem modelu referencyjnego h = const
Tb_ref = rozw_const(:,1);
err_Tb_poly = oblicz_bledy(Tb_ref, rozw_poly(:,1));
err_Tb_mnk  = oblicz_bledy(Tb_ref, rozw_mnk(:,1));
err_Tb_lin  = oblicz_bledy(Tb_ref, rozw_lin(:,1));

%% Wizualizacje
rysuj_h(DeltaT_pom, h_pom, DeltaT_siatka, h_poly, h_mnk_fit, h_lin);
rysuj_temperatury(t_const, rozw_const, t_poly, rozw_poly, t_mnk, rozw_mnk, t_lin, rozw_lin, param.h_const);
rysuj_bledy_Tb(t_const, rozw_poly(:,1), rozw_mnk(:,1), rozw_lin(:,1), Tb_ref);

%% Raport tekstowy
fprintf('\nTabela błędów h(DeltaT) względem pomiarów:\n');
drukuj_tabele_bledow({'Interpolacja poly','MNK','Interpolacja liniowa'}, [err_h_poly; err_h_mnk; err_h_lin]);

fprintf('\nTabela błędów Tb(t) względem modelu h=const:\n');
drukuj_tabele_bledow({'Interpolacja poly','MNK','Interpolacja liniowa'}, [err_Tb_poly; err_Tb_mnk; err_Tb_lin]);

fprintf('\nWnioski:\n');
interpretuj_wyniki(err_h_poly, err_h_mnk, err_h_lin, err_Tb_poly, err_Tb_mnk, err_Tb_lin);

end

%% --- Modele h(DeltaT) -------------------------------------------------
function h_val = h_interp_poly(DeltaT, DeltaT_pom, h_pom)
% interpolacja wielomianowa metodą Lagrange'a (z zestawu 1 zad. 4.2)
DeltaT = DeltaT(:).';
h_val = zeros(size(DeltaT));
for i = 1:numel(DeltaT)
    h_val(i) = lagrange(DeltaT_pom, h_pom, DeltaT(i));
end
h_val = h_val(:);
end

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

function h_val = h_interp_liniowa(DeltaT, DeltaT_pom, h_pom)
DeltaT = DeltaT(:).';
h_val = zeros(size(DeltaT));
for i = 1:numel(DeltaT)
    dt = DeltaT(i);
    if dt <= DeltaT_pom(1)
        h_val(i) = h_pom(1);
    elseif dt >= DeltaT_pom(end)
        h_val(i) = h_pom(end);
    else
        idx = find(DeltaT_pom <= dt, 1, 'last');
        if DeltaT_pom(idx) == dt
            h_val(i) = h_pom(idx);
        else
            dt1 = DeltaT_pom(idx); dt2 = DeltaT_pom(idx+1);
            h1 = h_pom(idx);      h2 = h_pom(idx+1);
            w = (dt - dt1) / (dt2 - dt1);
            h_val(i) = h1 + (h2 - h1) * w;
        end
    end
end
h_val = h_val(:);
end

%% --- Funkcje ODE ------------------------------------------------------
function [t, rozw] = symuluj_model(h_fun, param)
ode = @(czas, x, par) f_zmienny_h(czas, x, par, h_fun);
[t, rozw] = metodaEuleraUlepszona(ode, param.tspan(1), param.tspan(2), param.krok, [param.Tb0; param.Tw0], param);
end

function dx = f_zmienny_h(~, x, param, h_fun)
Tb = x(1); Tw = x(2);
h = h_fun(Tb, Tw);
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

%% --- Narzędzia matematyczne ------------------------------------------
function [err_struct, wsp, k_opt] = dopasuj_mnk(DeltaT, h)
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

function wart = wartosc_wielomianu(wsp, x)
wart = 0;
for i = 1:length(wsp)
    wart = wart * x + wsp(i);
end
end

function bledy = oblicz_bledy(ref, approx)
roznica = approx - ref;
bledy.MSE = mean(roznica.^2);
bledy.max = max(abs(roznica));
bledy.rel_mean = mean(abs(roznica)) / mean(abs(ref));
end

%% --- Wizualizacje i raporty ------------------------------------------
function rysuj_h(DeltaT_pom, h_pom, DeltaT_siatka, h_poly, h_mnk_fit, h_lin)
figure; hold on;
scatter(DeltaT_pom, h_pom, 40, 'k', 'filled', 'DisplayName','pomiary');
plot(DeltaT_siatka, h_poly, 'LineWidth',1.4,'DisplayName','interpolacja poly');
plot(DeltaT_siatka, h_mnk_fit,'LineWidth',1.4,'DisplayName','MNK');
plot(DeltaT_siatka, h_lin, 'LineWidth',1.4,'DisplayName','interpolacja liniowa');
plot(DeltaT_siatka, 160*ones(size(DeltaT_siatka)), 'k--','DisplayName','h const=160');
xlabel('\DeltaT [C]'); ylabel('h [W/(m^2*K)]');
title('Modele h(\DeltaT)'); grid on; legend('Location','best');
end

function rysuj_temperatury(t_const, rozw_const, t_poly, rozw_poly, t_mnk, rozw_mnk, t_lin, rozw_lin, h_const)
figure;
subplot(2,1,1); hold on;
plot(t_const, rozw_const(:,1),'k--','LineWidth',1.2,'DisplayName','h const');
plot(t_poly,  rozw_poly(:,1), 'LineWidth',1.2,'DisplayName','interpolacja poly');
plot(t_mnk,   rozw_mnk(:,1),  'LineWidth',1.2,'DisplayName','MNK');
plot(t_lin,   rozw_lin(:,1),  'LineWidth',1.2,'DisplayName','interpolacja liniowa');
title(sprintf('T_b(t) dla różnych h(\DeltaT) (h const = %g)', h_const));
xlabel('t [s]'); ylabel('T_b [C]'); grid on; legend('Location','best');

subplot(2,1,2); hold on;
plot(t_const, rozw_const(:,2),'k--','LineWidth',1.2,'DisplayName','h const');
plot(t_poly,  rozw_poly(:,2), 'LineWidth',1.2,'DisplayName','interpolacja poly');
plot(t_mnk,   rozw_mnk(:,2),  'LineWidth',1.2,'DisplayName','MNK');
plot(t_lin,   rozw_lin(:,2),  'LineWidth',1.2,'DisplayName','interpolacja liniowa');
title('T_w(t) dla różnych modeli h');
xlabel('t [s]'); ylabel('T_w [C]'); grid on; legend('Location','best');
end

function rysuj_bledy_Tb(t_ref, Tb_poly, Tb_mnk, Tb_lin, Tb_ref)
figure; hold on;
plot(t_ref, Tb_poly - Tb_ref, 'LineWidth',1.2,'DisplayName','interpolacja poly');
plot(t_ref, Tb_mnk  - Tb_ref, 'LineWidth',1.2,'DisplayName','MNK');
plot(t_ref, Tb_lin  - Tb_ref, 'LineWidth',1.2,'DisplayName','interpolacja liniowa');
plot(t_ref, zeros(size(t_ref)), 'k--','DisplayName','odniesienie');
title('Błąd T_b(t) względem modelu h=const');
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

function interpretuj_wyniki(err_h_poly, err_h_mnk, err_h_lin, err_Tb_poly, err_Tb_mnk, err_Tb_lin)
err_h = [err_h_poly.MSE, err_h_mnk.MSE, err_h_lin.MSE];
err_T = [err_Tb_poly.MSE, err_Tb_mnk.MSE, err_Tb_lin.MSE];
[~, idx_h] = min(err_h);
[~, idx_T] = min(err_T);
model_names = {'interpolacja wielomianowa','MNK','interpolacja liniowa'};
fprintf('Najlepsze dopasowanie h(DeltaT): %s\n', model_names{idx_h});
fprintf('Najmniejszy błąd przebiegu T_b(t): %s\n', model_names{idx_T});
if idx_h == idx_T
    fprintf('Spójny wybór modelu: %s jest rekomendowany do symulacji.\n', model_names{idx_h});
else
    fprintf('Różnice w kryteriach wskazują na potrzebę dodatkowej analizy (h: %s, T_b: %s).\n', model_names{idx_h}, model_names{idx_T});
end
end

%% --- Dane pomocnicze --------------------------------------------------
function [DeltaT_pom, h_pom] = pobierz_dane_pomiarowe()
DeltaT_pom = [-1500 -1000 -300 -50 -1 1 20 50 200 400 1000 2000];
h_pom      = [178   176   168  161 160 160 160.2 161 165 168 174 179];
end
