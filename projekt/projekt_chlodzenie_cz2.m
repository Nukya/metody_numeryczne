function projekt_chlodzenie_cz2
clc; clear; close all;

%% Parametry modelu fizycznego
parametry.h0 = 80;       % wartość bazowa współczynnika wnikania
parametry.alpha = 0.18;  % czułość na różnicę temperatur
parametry.beta = 0.25;   % nieliniowy wykładnik
parametry.A  = 0.0109;   % powierzchnia wymiany [m^2]
parametry.mb = 0.2;      % masa pręta [kg]
parametry.mw = 2.5;      % masa oleju [kg]
parametry.cb = 3.85;     % ciepło właściwe pręta [kJ/(kg*K)]
parametry.cw = 4.1813;   % ciepło właściwe oleju [kJ/(kg*K)]

%% Warunki początkowe
Tb0 = 1200;   % temperatura początkowa pręta [C]
Tw0 = 25;     % temperatura początkowa oleju [C]
stanPoczatkowy = [Tb0; Tw0];

%% Parametry numeryczne
przedzialCzasu = [0 5];
krokCzasu = 0.0025;               % krok czasu

%% 1. Symulacja – metoda Eulera
[czasE, rozwE] = metodaEulera(@f, przedzialCzasu(1), przedzialCzasu(2), krokCzasu, stanPoczatkowy, parametry);

%% 2. Symulacja – metoda zmodyfikowanego Eulera (Heun)
[czasM, rozwM] = metodaEuleraUlepszona(@f, przedzialCzasu(1), przedzialCzasu(2), krokCzasu, stanPoczatkowy, parametry);

%% 3. Symulacja ode45 – rozwiązanie referencyjne
odefun = @(t,x) f(t,x,parametry);
[czas45, rozw45] = ode45(odefun, przedzialCzasu, stanPoczatkowy);

%% 4. Wykresy Tb(t), Tw(t)
figure;
subplot(2,1,1); hold on;
plot(czasE, rozwE(:,1),'b','LineWidth',1.3,'DisplayName','Euler');
plot(czasM, rozwM(:,1),'r','LineWidth',1.3,'DisplayName','Mod. Euler');
plot(czas45, rozw45(:,1),'k--','LineWidth',1.3,'DisplayName','ode45');
title('Temperatura pręta T_b(t) – h nieliniowe');
xlabel('t [s]'); ylabel('T_b [C]'); legend; grid on;

subplot(2,1,2); hold on;
plot(czasE, rozwE(:,2),'b','LineWidth',1.3,'DisplayName','Euler');
plot(czasM, rozwM(:,2),'r','LineWidth',1.3,'DisplayName','Mod. Euler');
plot(czas45, rozw45(:,2),'k--','LineWidth',1.3,'DisplayName','ode45');
title('Temperatura oleju T_w(t) – h nieliniowe');
xlabel('t [s]'); ylabel('T_w [C]'); legend; grid on;

%% 5. Analiza wrażliwości na krok h (porównanie nieliniowe)
krokiTestowe = [0.01, 0.0025];

figure; hold on;
for krokTestu = krokiTestowe
    [czasTestu, rozwTestu] = metodaEulera(@f, przedzialCzasu(1), przedzialCzasu(2), krokTestu, stanPoczatkowy, parametry);
    plot(czasTestu, rozwTestu(:,1),'LineWidth',1.2,'DisplayName',['Euler h=',num2str(krokTestu)]);
end
plot(czas45, rozw45(:,1),'k--','LineWidth',1.3,'DisplayName','ode45');
title('Wpływ kroku h na rozwiązanie (h nieliniowe)'); xlabel('t'); ylabel('T_b');
legend; grid on;

%% 6. Porównanie: h stałe vs h nieliniowe
parametry_liniowe = parametry;
parametry_liniowe.alpha = 0; % brak zależności od temperatury
parametry_liniowe.beta  = 1;
parametry_liniowe.h0    = 160; % wartość z części 1

[czasStale, rozwStale] = metodaEuleraUlepszona(@f, przedzialCzasu(1), przedzialCzasu(2), krokCzasu, stanPoczatkowy, parametry_liniowe);

figure; hold on;
plot(czasM, rozwM(:,1),'r','LineWidth',1.3,'DisplayName','h nieliniowe');
plot(czasStale, rozwStale(:,1),'b--','LineWidth',1.3,'DisplayName','h stałe');
plot(czas45, rozw45(:,1),'k:','LineWidth',1.3,'DisplayName','ode45 nieliniowe');
title('Wpływ modelu h(T_b,T_w) na T_b(t)'); xlabel('t [s]'); ylabel('T_b [C]'); legend; grid on;

%% 7. Mapa h(Tb,Tw) – wizualizacja modelu
Tb_siatka = linspace(50, 1200, 60);
Tw_siatka = linspace(20, 200, 40);
[TB, TW] = meshgrid(Tb_siatka, Tw_siatka);
hMapa = fh(TB, TW, parametry);

figure;
surf(TB, TW, hMapa, 'EdgeColor','none');
colorbar;
title('Mapa współczynnika h(T_b,T_w)');
xlabel('T_b [C]'); ylabel('T_w [C]'); zlabel('h [W/(m^2*K)]');
view(45,30);

end

%% Funkcje pomocnicze

%% Równania stanu (ODE)
function dx = f(~, x, parametry)
Tb = x(1);
Tw = x(2);

wspWnikania = fh(Tb, Tw, parametry);
powierzchnia = parametry.A;
masaPreta = parametry.mb;
masaOleju = parametry.mw;
cieploPreta = parametry.cb;
cieploOleju = parametry.cw;

dTb = -(wspWnikania*powierzchnia)/(masaPreta*cieploPreta) * (Tb - Tw);
dTw =  +(wspWnikania*powierzchnia)/(masaOleju*cieploOleju) * (Tb - Tw);

dx = [dTb; dTw];
end

%% Model współczynnika h(Tb,Tw) – nieliniowy
function wsp_h = fh(Tb, Tw, parametry)
roznica = max(Tb - Tw, 0);
wsp_h = parametry.h0 + parametry.alpha * roznica .^ parametry.beta;
end
