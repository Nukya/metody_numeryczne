function projekt_chlodzenie_cz3
    clc; clear; close all;

    %%  1. Dane pomiarowe h(ΔT) 
    [dT_tab, h_tab] = dane_h();

    % Równoodległe węzły + h w tych punktach
    N = numel(dT_tab);
    dT_eq = linspace(min(dT_tab), max(dT_tab), N);
    h_eq  = wezly_rownoodlegle(dT_tab, h_tab, dT_eq);

    %%  2. Trzy modele h(ΔT) 

    % MNK
    stopien = 5;
    [wsp_mnk, ~] = aproksymacja_wielomianowa(dT_tab, h_tab, stopien);
    fh_mnk = @(dT) wartosc_wielomianu(wsp_mnk, dT);

    % Lagrange
    fh_lagrange = @(dT) lagrange(dT_eq, h_eq, dT);

    % Splajn – model referencyjny
    alpha = 0; beta = 0;
    S_spline = funkcje_sklejane(dT_eq, h_eq, alpha, beta);
    fh_spline = @(dT) S_spline.eval(dT);

    %%  3. Parametry materiałowe 
    parametry.A  = 0.0109;
    parametry.mb = 0.2;
    parametry.cb = 3.85;
    parametry.cw = 4.1813;

    t0 = 0;
    h = 0.001;

    %%  4. Tabela pomiarowa - wybór przypadków 
    dane = pomiary();
    
    % Wybieramy reprezentatywne przypadki: 1, 3, 7, 10
    idx = [1 3 7 10];
    P = dane(idx,:);
    n_przypadkow = length(idx);

    %%  5. Symulacje i analiza 
    
    % Tablice na wyniki
    bledy_Tb = zeros(n_przypadkow, 3); % [MNK, Lagrange, Splajn]
    bledy_Tw = zeros(n_przypadkow, 3);
    bledy_Tb_wzgl = zeros(n_przypadkow, 3);
    bledy_Tw_wzgl = zeros(n_przypadkow, 3);

    fprintf('CZĘŚĆ 3: ANALIZA METOD APROKSYMACJI h(ΔT)\n');

    for k = 1:n_przypadkow
        Tb0 = P(k,1);
        Tw0 = P(k,2);
        parametry.mw = P(k,3);
        tK  = P(k,4);
        x0  = [Tb0; Tw0];

        % Symulacje dla trzech metod
        [tM, xM] = metodaEuleraUlepszona_nielin(@f_nielin, t0, tK, h, x0, parametry, fh_mnk);
        [tL, xL] = metodaEuleraUlepszona_nielin(@f_nielin, t0, tK, h, x0, parametry, fh_lagrange);
        [tS, xS] = metodaEuleraUlepszona_nielin(@f_nielin, t0, tK, h, x0, parametry, fh_spline);

        %% Wartości końcowe
        Tb_meas = P(k,5);
        Tw_meas = P(k,6);

        Tb_MNK = xM(end,1);
        Tb_LAG = xL(end,1);
        Tb_SPL = xS(end,1);

        Tw_MNK = xM(end,2);
        Tw_LAG = xL(end,2);
        Tw_SPL = xS(end,2);

        %% Obliczenie błędów
        bledy_Tb(k,:) = [Tb_MNK - Tb_meas, Tb_LAG - Tb_meas, Tb_SPL - Tb_meas];
        bledy_Tw(k,:) = [Tw_MNK - Tw_meas, Tw_LAG - Tw_meas, Tw_SPL - Tw_meas];
        
        bledy_Tb_wzgl(k,:) = abs(bledy_Tb(k,:)) ./ Tb_meas * 100;
        bledy_Tw_wzgl(k,:) = abs(bledy_Tw(k,:)) ./ Tw_meas * 100;

        %% Wyświetlenie wyników
        fprintf('\n===== Zestaw danych %d (Tb0=%.0f, Tw0=%.0f, mw=%.1f, t=%.0fs) =====\n', ...
            idx(k), Tb0, Tw0, parametry.mw, tK);
        fprintf('Pomiar:   Tb = %.2f°C,  Tw = %.2f°C\n', Tb_meas, Tw_meas);
        fprintf('MNK:      Tb = %.2f°C,  Tw = %.2f°C  (ΔTb = %+.2f, ΔTw = %+.2f)\n', ...
            Tb_MNK, Tw_MNK, bledy_Tb(k,1), bledy_Tw(k,1));
        fprintf('Lagrange: Tb = %.2f°C,  Tw = %.2f°C  (ΔTb = %+.2f, ΔTw = %+.2f)\n', ...
            Tb_LAG, Tw_LAG, bledy_Tb(k,2), bledy_Tw(k,2));
        fprintf('Splajn:   Tb = %.2f°C,  Tw = %.2f°C  (ΔTb = %+.2f, ΔTw = %+.2f)\n', ...
            Tb_SPL, Tw_SPL, bledy_Tb(k,3), bledy_Tw(k,3));
        
        %% Wykresy dla każdego przypadku
        figure('Name', sprintf('porownanie metod - przypadek %d', idx(k)));
        
        % Tb(t)
        subplot(2,1,1);
        hold on; grid on;
        plot(tM, xM(:,1), 'Color',[1 .5 0],'LineWidth',1.5, 'DisplayName','MNK');
        plot(tL, xL(:,1), 'r','LineWidth',1.5, 'DisplayName','Lagrange');
        plot(tS, xS(:,1), 'm','LineWidth',1.5, 'DisplayName','Splajn');
        scatter(tK, Tb_meas, 100, 'k', 'filled', 'DisplayName','Pomiar');
        title(sprintf('T_b(t) – Przypadek %d', idx(k)));
        xlabel('t [s]'); ylabel('T_b [°C]');
        legend('Location','best');
        
        % Tw(t)
        subplot(2,1,2);
        hold on; grid on;
        plot(tM, xM(:,2), 'Color',[1 .5 0],'LineWidth',1.5, 'DisplayName','MNK');
        plot(tL, xL(:,2), 'r','LineWidth',1.5, 'DisplayName','Lagrange');
        plot(tS, xS(:,2), 'm','LineWidth',1.5, 'DisplayName','Splajn');
        scatter(tK, Tw_meas, 100, 'k','filled', 'DisplayName','Pomiar');
        title(sprintf('T_w(t) – Przypadek %d', idx(k)));
        xlabel('t [s]'); ylabel('T_w [°C]');
        legend('Location','best');
    end

    %% 6. ANALIZA PORÓWNAWCZA - Podsumowanie
    fprintf('PODSUMOWANIE BŁĘDÓW\n');
    
    % Błędy bezwzględne średnie
    bledy_Tb_mean = mean(abs(bledy_Tb), 1);
    bledy_Tw_mean = mean(abs(bledy_Tw), 1);
    
    % Błędy maksymalne
    bledy_Tb_max = max(abs(bledy_Tb), [], 1);
    bledy_Tw_max = max(abs(bledy_Tw), [], 1);
    
    % Błędy względne średnie
    bledy_Tb_wzgl_mean = mean(bledy_Tb_wzgl, 1);
    bledy_Tw_wzgl_mean = mean(bledy_Tw_wzgl, 1);
    
    fprintf('\nŚrednie błędy bezwzględne temperatury pręta Tb:\n');
    fprintf('  MNK:      %.2f°C\n', bledy_Tb_mean(1));
    fprintf('  Lagrange: %.2f°C\n', bledy_Tb_mean(2));
    fprintf('  Splajn:   %.2f°C\n', bledy_Tb_mean(3));
    
    fprintf('\nMaksymalne błędy bezwzględne temperatury pręta Tb:\n');
    fprintf('  MNK:      %.2f°C\n', bledy_Tb_max(1));
    fprintf('  Lagrange: %.2f°C\n', bledy_Tb_max(2));
    fprintf('  Splajn:   %.2f°C\n', bledy_Tb_max(3));
    
    fprintf('\nŚrednie błędy względne temperatury pręta Tb:\n');
    fprintf('  MNK:      %.2f%%\n', bledy_Tb_wzgl_mean(1));
    fprintf('  Lagrange: %.2f%%\n', bledy_Tb_wzgl_mean(2));
    fprintf('  Splajn:   %.2f%%\n', bledy_Tb_wzgl_mean(3));
    
    fprintf('\nŚrednie błędy bezwzględne temperatury oleju Tw:\n');
    fprintf('  MNK:      %.2f°C\n', bledy_Tw_mean(1));
    fprintf('  Lagrange: %.2f°C\n', bledy_Tw_mean(2));
    fprintf('  Splajn:   %.2f°C\n', bledy_Tw_mean(3));
    
    fprintf('\nŚrednie błędy względne temperatury oleju Tw:\n');
    fprintf('  MNK:      %.2f%%\n', bledy_Tw_wzgl_mean(1));
    fprintf('  Lagrange: %.2f%%\n', bledy_Tw_wzgl_mean(2));
    fprintf('  Splajn:   %.2f%%\n', bledy_Tw_wzgl_mean(3));
    
end