function projekt_chlodzenie_cz2
    clc; clear; close all;

    
    %  1. Dane pomiarowe h(ΔT)
    
    [dT_tab, h_tab] = dane_h();


        % === Wygenerowanie równoodległych węzłów ΔT ===
    N = numel(dT_tab);
    dT_eq = linspace(min(dT_tab), max(dT_tab), N);
    
    %  Metoda 1: aproksymacja MNK 
    stopien = 5;
    [wsp_mnk, ~] = aproksymacja_wielomianowa(dT_tab, h_tab, stopien);
    h_eq = interpolacja_liniowa(dT_tab, h_tab, dT_eq);    
    fh_mnk = @(dT) wartosc_wielomianu(wsp_mnk, dT);

    %  2. Metoda 2 – interpolacja Lagrange'a
    fh_lagrange = @(dT) lagrange(dT_eq, h_eq, dT);
    
    %  4. Metoda 3 – funkcje sklejane (splajn kubiczny)

    % warunki brzegowe (pochodne)
    alpha = 0;
    beta  = 0;

    S_spline = funkcje_sklejane(dT_eq, h_eq, alpha, beta);
    fh_spline = @(dT) S_spline.eval(dT);



    %  5. Wykres porównawczy trzech metod h(ΔT)

    dT_dense = linspace(min(dT_tab), max(dT_tab), 500);
         
    %% --- 1. MNK ---
    figure;
    hold on; grid on;
    plot(dT_tab, h_tab, 'ko', 'MarkerFaceColor','k', 'DisplayName','Dane pomiarowe');
    plot(dT_dense, fh_mnk(dT_dense), 'Color',[1 0.5 0], 'LineWidth',1.5, 'DisplayName','MNK');
    
    title('Aproksymacja MNK');
    xlabel('\Delta T [°C]');
    ylabel('h [W/m^2]');
    legend('Location','best');

    %% --- 2. Lagrange ---
    figure;
    hold on; grid on;
    plot(dT_tab, h_tab, 'ko', 'MarkerFaceColor','k', 'DisplayName','Dane pomiarowe');
    plot(dT_dense, fh_lagrange(dT_dense), 'r', 'LineWidth',1.5, 'DisplayName','Lagrange');
    
    title('Interpolacja Lagrange''a');
    xlabel('\Delta T [°C]');
    ylabel('h [W/m^2]');
    legend('Location','best');
    
    
    
    %% --- 3. Splajn ---
    figure;
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
    [tL,  xL]  = metodaEuleraUlepszona_nielin(@f_nielin, t0, tK, h, x0, parametry, fh_lagrange);
    [tS,  xS]  = metodaEuleraUlepszona_nielin(@f_nielin, t0, tK, h, x0, parametry, fh_spline);


    %  7. Wykres TB(t) i TW(t)

    figure;

    hold on; 
    grid on;
    plot(tM, xM(:,1), 'LineWidth',1.3, 'DisplayName','MNK');
    plot(tS, xS(:,1), 'LineWidth',1.3, 'DisplayName','Splajn');
    plot(tL, xL(:,1), 'LineWidth',1.3, 'DisplayName','Lagrange');

        hold off;
        xlabel('t [s]'); 
        ylabel('T_b [°C]');
        title('T_b(t) dla omawianych metod');
        legend('Location','best');

    figure;
    hold on; 
    grid on;
    plot(tL, xL(:,2), 'LineWidth',1.3, 'DisplayName','MNK');
    plot(tS, xS(:,2), 'LineWidth',1.3, 'DisplayName','Splajn');
    plot(tM, xM(:,2), 'LineWidth',1.3, 'DisplayName','Lagrange');

        hold off;
        xlabel('t [s]'); 
        ylabel('T_w [°C]');
        title('T_w(t) dla omawianych metod');
        legend('Location','best');

    %% Obliczenie błędu bezwzględnego i względnego
    h_mnk_tab      = fh_mnk(dT_tab);
    h_lagrange_tab = fh_lagrange(dT_tab);
    h_spline_tab   = fh_spline(dT_tab);
    
    metody = {'MNK','Lagrange','Splajn'};
    H = [h_mnk_tab(:), h_lagrange_tab(:), h_spline_tab(:)];
    
    blad_bezwzgledny_sredni = zeros(3,1);
    blad_bezwzgledny_max = zeros(3,1);
    
    blad_wzgledny_sredni = zeros(3,1);
    blad_wzgledny_max = zeros(3,1);
    
    for i = 1:3
        err_abs = abs(H(:,i) - h_tab(:));                    % błąd bezwzględny
        err_rel = abs((H(:,i) - h_tab(:)) ./ h_tab(:)) * 100; % błąd względny [%]
    
        blad_bezwzgledny_sredni(i) = mean(err_abs);
        blad_bezwzgledny_max(i)    = max(err_abs);
    
        blad_wzgledny_sredni(i) = mean(err_rel);
        blad_wzgledny_max(i)    = max(err_rel);
    end
    
    T_bledy = table( ...
        metody', ...
        blad_bezwzgledny_sredni, ...
        blad_bezwzgledny_max, ...
        blad_wzgledny_sredni, ...
        blad_wzgledny_max, ...
        'VariableNames', { ...
            'Metoda', ...
            'blad_bezwzgledny_sredni', ...
            'blad_bezwzgledny_max', ...
            'blad_wzgledny_sredni', ...
            'blad_wzgledny_max' ...
        } ...
    );
    
    disp(' ');
    disp('===== Błędy bezwzględne i względne =====');
    disp(T_bledy);



end
