function projekt_chlodzenie_cz1
    clc; 
    clear; 
    close all;
    
    %% Parametry modelu fizycznego
    parametry.h  = 160;      % wsp. wnikania ciepła [W/(m^2*K)]
    parametry.A  = 0.0109;   % powierzchnia wymiany [m^2]
    parametry.mb = 0.2;      % masa pręta [kg]
    parametry.mw = 2.5;      % masa oleju [kg]
    parametry.cb = 3.85;     % ciepło właśc. pręta [kJ/(kg*K)]
    parametry.cw = 4.1813;   % ciepło właśc. oleju [kJ/(kg*K)]
    
    
    %% Warunki początkowe (konfiguracja bazowa)
    Tb0 = 1200;   % temperatura początkowa pręta [C]
    Tw0 = 25;     % temperatura początkowa oleju [C]
    stanPoczatkowy = [Tb0; Tw0];
    
    %% Parametry numeryczne
    przedzialCzasu = [0 5];
    krokCzasu      = 0.001;   % krok czasu (do symulacji bazowej)
    
    %% 1. Symulacja – metoda Eulera
    [czasE, rozwE] = metodaEulera(@f, przedzialCzasu(1), przedzialCzasu(2), ...
                                  krokCzasu, stanPoczatkowy, parametry);
    
    %% 2. Symulacja – metoda zmodyfikowanego Eulera (Heun)
    [czasM, rozwM] = metodaEuleraUlepszona(@f, przedzialCzasu(1), przedzialCzasu(2), ...
                                           krokCzasu, stanPoczatkowy, parametry);
    
    %% 3. Symulacja ode45 – rozwiązanie referencyjne
    odefun = @(t,x) f(t,x,parametry);
    [czas45, rozw45] = ode45(odefun, przedzialCzasu, stanPoczatkowy);
    
    %% 4. Wykresy Tb(t), Tw(t) – porównanie metod
    figure('Name', 'przebieg temperatur dla roznych metod rozwiazywania rownan rozniczkowych');
    subplot(2,1,1); 
    hold on; 
    grid on;
    plot(czasE,  rozwE(:,1),'b','LineWidth',1.3,'DisplayName','Euler');
    plot(czasM,  rozwM(:,1),'r','LineWidth',1.3,'DisplayName','Mod. Euler');
    plot(czas45, rozw45(:,1),'k--','LineWidth',1.3,'DisplayName','ode45');
        title('Temperatura pręta T_b(t)');
        xlabel('t [s]'); 
        ylabel('T_b [C]'); 
        legend('Location','best');
    
    subplot(2,1,2); 
    hold on; 
    grid on;
    plot(czasE,  rozwE(:,2),'b','LineWidth',1.3,'DisplayName','Euler');
    plot(czasM,  rozwM(:,2),'r','LineWidth',1.3,'DisplayName','Mod. Euler');
    plot(czas45, rozw45(:,2),'k--','LineWidth',1.3,'DisplayName','ode45');
        title('Temperatura oleju T_w(t)');
        xlabel('t [s]'); 
        ylabel('T_w [C]'); 
        legend('Location','best');
    
    %% 5. Analiza wrażliwości na krok h (EULER vs ode45)
    krokiTestowe = [0.2, 0.1, 0.001];
    
    figure('Name', 'wplyw kroku na rozwiazanie - euler'); 
    hold on; 
    grid on;
    for krokTestu = krokiTestowe
        [czasTestu, rozwTestu] = metodaEulera(@f, przedzialCzasu(1), przedzialCzasu(2), ...
                                              krokTestu, stanPoczatkowy, parametry);
        plot(czasTestu, rozwTestu(:,1),'LineWidth',1.2, ...
             'DisplayName',['Euler, h = ',num2str(krokTestu)]);
    end
    plot(czas45, rozw45(:,1),'k--','LineWidth',1.3,'DisplayName','ode45');
        title('Wpływ kroku h na rozwiązanie - metoda Eulera');
        xlabel('t [s]'); 
        ylabel('T_b [C]');
        legend('Location','best');
    
    %% 6. Analiza wrażliwości na krok h (ULEPSZONY EULER vs ode45)
    figure('Name', 'wplyw kroku na rozwiazanie - euler ulepszony'); 
    hold on; 
    grid on;
    for krokTestu = krokiTestowe
        [czasTestu, rozwTestu] = metodaEuleraUlepszona(@f, przedzialCzasu(1), przedzialCzasu(2), ...
                                                       krokTestu, stanPoczatkowy, parametry);
        plot(czasTestu, rozwTestu(:,1),'LineWidth',1.2, ...
             'DisplayName',['Mod. Euler, h = ',num2str(krokTestu)]);
    end
    
    plot(czas45, rozw45(:,1),'k--','LineWidth',1.3,'DisplayName','ode45');
        title('Wpływ kroku h na rozwiązanie - metoda ulepszona Eulera');
        xlabel('t [s]'); 
        ylabel('T_b [C]');
        legend('Location','best');
    
    %% 7. Weryfikacja z danymi pomiarowymi + przebiegi dla wszystkich przypadków
    
    % Dane:         Tb0 Tw0 mw t_obs Tb_exp Tw_exp
    pomiary =   [
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
    
    %% 7. Weryfikacja z danymi pomiarowymi – każdy zestaw w osobnym wykresie
    Tb_model = zeros(10,1);
    Tw_model = zeros(10,1);
    dane = zeros(10,7);   % [Tb_mod Tb_exp Tw_mod Tw_exp ErrTb ErrTw]
    
    figure('Name', 'przebiegi dla roznych zestawow danych');
    t = tiledlayout(5,2);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';


    for i = 1:10
        parametry_i = parametry;
        parametry_i.mw = pomiary(i,3);      % masa oleju z tabeli
    
        Tb0_i = pomiary(i,1);               %temperatura początkowa pręta
        Tw0_i = pomiary(i,2);               %temperatura początkowa cieczy chłodzącej
        t_obs = pomiary(i,4);               %czas obserwacji
        Tb_exp = pomiary(i,5);
        Tw_exp = pomiary(i,6);
    
        [czasLokalny, rozwLokalne] = metodaEuleraUlepszona(@f, 0, t_obs, ...
            krokCzasu, [Tb0_i; Tw0_i], parametry_i);
    
        Tb_model(i) = rozwLokalne(end,1);
        Tw_model(i) = rozwLokalne(end,2);
    
        dane(i,:) = [i ,Tb_model(i), Tb_exp, Tw_model(i), Tw_exp, ...
                     Tb_model(i)-Tb_exp, Tw_model(i)-Tw_exp];
    
        nexttile;
        hold on;
        grid on;
        plot(czasLokalny, rozwLokalne(:,1), 'r', 'LineWidth', 1.3, ...
            'DisplayName', 'T_b(t)');
        plot(czasLokalny, rozwLokalne(:,2), 'b', 'LineWidth', 1.3, ...
            'DisplayName', 'T_w(t)');
    
            title(sprintf('Zestaw danych %d', i));
            xlabel('t [s]');
            ylabel('Temperatura [°C]');
            legend('Location','best');
    
    end
    
    fprintf("\n WERYFIKACJA MODEL–POMIAR \n");
    T = array2table(dane, ...
        'VariableNames', {'nr', 'Tb_mod','Tb_exp','Tw_mod','Tw_exp','ErrTb','ErrTw'});
    
    disp(T);
    
    errTb = T.ErrTb;
    errTw = T.ErrTw;
    [maxErrTb, idxTb] = max(abs(errTb));
    [maxErrTw, idxTw] = max(abs(errTw));
    
    fprintf("Największy błąd Tw = %.4f w zestawie %d\n", maxErrTw, idxTw);
    fprintf("Największy błąd Tb = %.4f w zestawie %d\n", maxErrTb, idxTb);
    
    
    %% 8. Wrażliwość na błędy danych wejściowych
    Tb0_bazowe = 1200;
    Tw0_bazowe = 25;
    mw_bazowe  = 2.5;
    czasObserwacji = 3;
    
    Tb_odniesienie = symuluj(Tb0_bazowe,Tw0_bazowe,mw_bazowe, ...
                                         czasObserwacji,parametry,krokCzasu);
    
    przypadki = {
        'Tb0+10', Tb0_bazowe+10, Tw0_bazowe,   mw_bazowe;
        'Tb0-10', Tb0_bazowe-10, Tw0_bazowe,   mw_bazowe;
        'Tw0+10',  Tb0_bazowe,    Tw0_bazowe+10, mw_bazowe;
        'Tw0-10',  Tb0_bazowe,    Tw0_bazowe-10, mw_bazowe;
        'mw+5%',  Tb0_bazowe,    Tw0_bazowe,   mw_bazowe*1.05;
        'mw-5%',  Tb0_bazowe,    Tw0_bazowe,   mw_bazowe*0.95;
    };
    
    liczbaPrzypadkow = size(przypadki,1);
    deltaTb = zeros(liczbaPrzypadkow,1);
    Tb_koniec_przypadki = zeros(liczbaPrzypadkow,1);
    
    for i = 1:liczbaPrzypadkow
        Tb_koniec_przypadki(i) = symuluj(przypadki{i,2}, przypadki{i,3}, przypadki{i,4}, ...
            czasObserwacji, parametry, krokCzasu);
        deltaTb(i) = Tb_koniec_przypadki(i) - Tb_odniesienie;
    end
    
    fprintf("\n WRAŻLIWOŚĆ NA BŁĘDY DANYCH WEJŚCIOWYCH \n");
    fprintf("\nReferencyjne T_b(t=%.1f s): %.4f C\n\n", czasObserwacji, Tb_odniesienie);
    
    etykieta = string(przypadki(:,1));
    wynikiWrazliwosci = table(etykieta, Tb_koniec_przypadki, deltaTb, ...
        'VariableNames', {'Błąd', 'Tb_koniec', 'DeltaTb'});
    
    disp(wynikiWrazliwosci);
    
end