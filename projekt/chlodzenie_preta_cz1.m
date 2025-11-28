function chlodzenie_preta_cz1
clc; clear; close all;

%% Parametry modelu fizycznego
parametry.h  = 160;
parametry.A  = 0.0109;
parametry.mb = 0.2;
parametry.mw = 2.5;
parametry.cb = 3.85;
parametry.cw = 4.1813;

%% Warunki początkowe
Tb0 = 1200;
Tw0 = 25;
stanPoczatkowy = [Tb0; Tw0];

%% Parametry numeryczne
przedzialCzasu = [0 5];
krokCzasu = 0.0025;        % wybrany „lepszy" krok

%   1. SYMULACJA – EULER
[czasE, rozwE] = metodaEulera(@f, przedzialCzasu(1), przedzialCzasu(2), krokCzasu, stanPoczatkowy, parametry);

%   2. SYMULACJA – ZMODYFIKOWANY EULER
[czasM, rozwM] = metodaEuleraUlepszona(@f, przedzialCzasu(1), przedzialCzasu(2), krokCzasu, stanPoczatkowy, parametry);

%   3. SYMULACJA – ODE45 (referencja)
odefun = @(t,x) f(t,x,parametry);
[czas45, rozw45] = ode45(odefun, przedzialCzasu, stanPoczatkowy);

%   4. WYKRESY: Tb(t), Tw(t)
figure;
subplot(2,1,1); hold on;
plot(czasE, rozwE(:,1),'b','LineWidth',1.3,'DisplayName','Euler');
plot(czasM, rozwM(:,1),'r','LineWidth',1.3,'DisplayName','Mod. Euler');
plot(czas45, rozw45(:,1),'k--','LineWidth',1.3,'DisplayName','ode45');
title('Temperatura pręta T_b(t)');
xlabel('t [s]'); ylabel('T_b [C]'); legend; grid on;

subplot(2,1,2); hold on;
plot(czasE, rozwE(:,2),'b','LineWidth',1.3,'DisplayName','Euler');
plot(czasM, rozwM(:,2),'r','LineWidth',1.3,'DisplayName','Mod. Euler');
plot(czas45, rozw45(:,2),'k--','LineWidth',1.3,'DisplayName','ode45');
title('Temperatura oleju T_w(t)');
xlabel('t [s]'); ylabel('T_w [C]'); legend; grid on;


%   5. ANALIZA WRAŻLIWOŚCI NA KROK h
krokiTestowe = [0.01, 0.0025];

figure; hold on;
for krokTestu = krokiTestowe
    [czasTestu, rozwTestu] = metodaEulera(@f, przedzialCzasu(1), przedzialCzasu(2), krokTestu, stanPoczatkowy, parametry);
    plot(czasTestu, rozwTestu(:,1),'LineWidth',1.2,'DisplayName',['Euler h=',num2str(krokTestu)]);
end
plot(czas45, rozw45(:,1),'k--','LineWidth',1.3,'DisplayName','ode45');
title('Wpływ kroku h na rozwiązanie'); xlabel('t'); ylabel('T_b');
legend; grid on;


%   6. WERYFIKACJA Z DANYMI POMIAROWYMI
pomiary = [
    1200 25 2.5 3 107.7 105.1
    800  25 2.5 3 79.1  78.0
    1100 70 2.5 3 142.1 139.1
    1200 25 2.5 5 105.7 105.5
    800  25 2.5 5 78.2  78.1
    1100 70 2.5 2 150.1 138.2
    1100 70 5   2 116.6 105.1
    1100 70 10  2 99.1  88.1
    1100 70 2.5 4 141.2 139.8
    1100 70 2.5 5 140.9 140.1
];

fprintf("Nr | Tb_mod | Tb_exp | Tw_mod | Tw_exp | ErrTb | ErrTw\n");

Tb_model = zeros(10,1);
Tw_model = zeros(10,1);

for i = 1:10

    parametry_i = parametry;
    parametry_i.mw = pomiary(i,3);   % zmiana masy oleju

    [czasLokalny, rozwLokalne] = metodaEuleraUlepszona(@f, 0, pomiary(i,4), krokCzasu, pomiary(i,1:2)', parametry_i);

    Tb_model(i) = rozwLokalne(end,1);
    Tw_model(i) = rozwLokalne(end,2);

    fprintf("%2d | %7.2f | %7.2f | %7.2f | %7.2f | %+6.2f | %+6.2f\n",...
        i, Tb_model(i), pomiary(i,5), Tw_model(i), pomiary(i,6), ...
        Tb_model(i)-pomiary(i,5), Tw_model(i)-pomiary(i,6));
end


%   7. WRAŻLIWOŚĆ NA BŁĘDY DANYCH WEJŚCIOWYCH

Tb0_bazowe = 1200;
Tw0_bazowe = 25;
mw_bazowe  = 2.5;
czasObserwacji = 3;

Tb_odniesienie = symuluj_jednorazowo(Tb0_bazowe,Tw0_bazowe,mw_bazowe,czasObserwacji,parametry,krokCzasu);

przypadki = {
    'Tb0+10', Tb0_bazowe+10, Tw0_bazowe,   mw_bazowe;
    'Tb0-10', Tb0_bazowe-10, Tw0_bazowe,   mw_bazowe;
    'Tw0+1',  Tb0_bazowe,    Tw0_bazowe+1, mw_bazowe;
    'Tw0-1',  Tb0_bazowe,    Tw0_bazowe-1, mw_bazowe;
    'mw+5%',  Tb0_bazowe,    Tw0_bazowe,   mw_bazowe*1.05;
    'mw-5%',  Tb0_bazowe,    Tw0_bazowe,   mw_bazowe*0.95;
};

liczbaPrzypadkow = size(przypadki,1);
deltaTb = zeros(liczbaPrzypadkow,1);

for i = 1:liczbaPrzypadkow
    deltaTb(i) = symuluj_jednorazowo(przypadki{i,2}, przypadki{i,3}, przypadki{i,4}, ...
        czasObserwacji, parametry, krokCzasu) - Tb_odniesienie;
end

figure;
bar(deltaTb); xticklabels(przypadki(:,1));
title('Wrażliwość na błędy danych wejściowych');
ylabel('\Delta T_b [C]');
grid on;

end

%% Funkcje pomocnicze

%% Równania stanu (ODE)
function dx = f(~, x, parametry)

Tb = x(1);
Tw = x(2);

wspWnikania = parametry.h;
powierzchnia = parametry.A;
masaPreta = parametry.mb;
masaOleju = parametry.mw;
cieploPreta = parametry.cb;
cieploOleju = parametry.cw;

dTb = -(wspWnikania*powierzchnia)/(masaPreta*cieploPreta) * (Tb - Tw);
dTw =  +(wspWnikania*powierzchnia)/(masaOleju*cieploOleju) * (Tb - Tw);

dx = [dTb; dTw];

end


%   MODEL WSPÓŁCZYNNIKA h(Tb,Tw)
%   (W CZĘŚCI 1 STAŁY – docelowo w Części 2 będzie nieliniowy)
function wsp_h = fh(Tb, Tw) %#ok<INUSD>
wsp_h = 160; % stały dla Części 1
end


%   SYMULACJA DO t_obs – używana w wrażliwości
function Tb_koniec = symuluj_jednorazowo(Tb0,Tw0,mw,czasObserwacji,parametry,krokCzasu)

parametryLokalne = parametry;
parametryLokalne.mw = mw;

[~, rozw] = metodaEuleraUlepszona(@f, 0, czasObserwacji, krokCzasu, [Tb0;Tw0], parametryLokalne);

Tb_koniec = rozw(end,1);

end
