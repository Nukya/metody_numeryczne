function projekt_chlodzenie_cz2
    clc; clear; close all;

    
    %  1. Dane pomiarowe h(ΔT)
    
    [dT_tab, h_tab] = dane_h();


    
    %  2. Metoda 1 – interpolacja Lagrange'a
    
    N = numel(dT_tab);
    
    % równoodległe węzły
    dT_eq = linspace(min(dT_tab), max(dT_tab), N);
    
    % wartości h w równoodległych węzłach
    h_eq = lagrange(dT_tab, h_tab, dT_eq);
    
    % interpolacja Lagrange'a na równoodległych węzłach
    fh_lagrange = @(dT) lagrange(dT_eq, h_eq, dT);


    
    %  3. Metoda 2 – aproksymacja MNK
    
    stopien = 5;
    [wsp_mnk,~] = aproksymacja_wielomianowa(dT_tab, h_tab, stopien);
    fh_mnk = @(dT) wartosc_wielomianu(wsp_mnk, dT);


    
    %  4. Metoda 3 – funkcje sklejane (splajn kubiczny)

    % warunki brzegowe (pochodne)
    alpha = 0;
    beta  = 0;

    S_spline = funkcje_sklejane(dT_eq, h_eq, alpha, beta);
    fh_spline = @(dT) S_spline.eval(dT);



    %  5. Wykres porównawczy trzech metod h(ΔT)

    dT_dense = linspace(min(dT_tab), max(dT_tab), 500);
    
    figure;
    tiledlayout(3,1);
    
    %% --- 1. Lagrange ---
    nexttile;
    hold on; grid on;
    plot(dT_tab, h_tab, 'ko', 'MarkerFaceColor','k', 'DisplayName','Dane pomiarowe');
    plot(dT_dense, fh_lagrange(dT_dense), 'r', 'LineWidth',1.5, 'DisplayName','Lagrange');
    
    title('Interpolacja Lagrange''a');
    xlabel('\Delta T [°C]');
    ylabel('h [W/m^2]');
    legend('Location','best');
    
    
    %% --- 2. MNK ---
    nexttile;
    hold on; grid on;
    plot(dT_tab, h_tab, 'ko', 'MarkerFaceColor','k', 'DisplayName','Dane pomiarowe');
    plot(dT_dense, fh_mnk(dT_dense), 'Color',[1 0.5 0], 'LineWidth',1.5, 'DisplayName','MNK');
    
    title('Aproksymacja MNK');
    xlabel('\Delta T [°C]');
    ylabel('h [W/m^2]');
    legend('Location','best');
    
    
    %% --- 3. Splajn ---
    nexttile;
    hold on; grid on;
    plot(dT_tab, h_tab, 'ko', 'MarkerFaceColor','k', 'DisplayName','Dane pomiarowe');
    plot(dT_dense, fh_spline(dT_dense), 'm', 'LineWidth',1.5, 'DisplayName','Splajn');
    
    title('Splajn kubiczny');
    xlabel('\Delta T [°C]');
    ylabel('h [W/m^2]');
    legend('Location','best');

    
    %  6. Symulacja chłodzenia
    
    parametry.A  = 0.0109;
    parametry.mb = 0.2;
    parametry.mw = 2.5;
    parametry.cb = 3.85;
    parametry.cw = 4.1813;

    Tb0 = 1200;
    Tw0 = 25;
    x0 = [Tb0; Tw0];

    t0 = 0; 
    tK = 5;
    h = 0.001;


    %% --- Symulacje dla wybranej metody ---
    [tM,  xM]  = metodaEuleraUlepszona_nielin(@f_nielin, t0, tK, h, x0, parametry, fh_mnk);

    %  7. Wykres TB(t) i TW(t)

    figure;

    hold on; 
    grid on;
    plot(tM, xM(:,1), 'LineWidth',1.3, 'DisplayName','T_b');
    plot(tM, xM(:,2), 'LineWidth',1.3, 'DisplayName','T_w');

        xlabel('t [s]'); 
        ylabel('T_w, T_b [°C]');
        title('T_w(t), t_b(t)');
        legend('Location','best');

end
