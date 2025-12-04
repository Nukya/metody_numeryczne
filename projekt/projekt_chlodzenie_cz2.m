function projekt_chlodzenie_cz2
    clc; clear; close all;

    %% 1. Dane pomiarowe h(ΔT)
    [dT_tab, h_tab] = dane_h();

    %% 2. Wygenerowanie równoodległych węzłów ΔT 
    N = numel(dT_tab);
    dT_eq = linspace(min(dT_tab), max(dT_tab), N);
    h_eq = wezly_rownoodlegle(dT_tab, h_tab, dT_eq);    

    %% 3.1 Metoda 1: Aproksymacja MNK 
    % MNK używa wszystkich oryginalnych punktów pomiarowych
    stopien = 5;
    [wsp_mnk, ~] = aproksymacja_wielomianowa(dT_tab, h_tab, stopien);
    fh_mnk = @(dT) wartosc_wielomianu(wsp_mnk, dT);

    %% 3.2 Metoda 2: Interpolacja Lagrange'a
    % Lagrange używa równoodległych węzłów (liczba węzłów = liczba danych)
    fh_lagrange = @(dT) lagrange(dT_eq, h_eq, dT);
    
    %% 3.3 Metoda 3: Funkcje sklejane (splajn kubiczny)
    % Splajn używa równoodległych węzłów (liczba węzłów = liczba danych)
    alpha = 0;  % Warunek brzegowy - pochodna na lewym końcu
    beta  = 0;  % Warunek brzegowy - pochodna na prawym końcu
    S_spline = funkcje_sklejane(dT_eq, h_eq, alpha, beta);
    fh_spline = @(dT) S_spline.eval(dT);
    
    fprintf('Liczba węzłów wynosi: %d\n', N);

    %% 6. Wykresy porównawcze trzech metod h(ΔT)
    % Gęsta siatka tylko do wizualizacji krzywych
    dT_dense = linspace(min(dT_tab), max(dT_tab), 500);
         
    % 6.1 MNK 
    figure('Name', 'Aproksymacja MNK');
    hold on; grid on;
    plot(dT_tab, h_tab, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Dane pomiarowe');
    plot(dT_dense, fh_mnk(dT_dense), 'Color', [1 0.5 0], 'LineWidth', 2, 'DisplayName', 'MNK');
    title('Aproksymacja MNK', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('\DeltaT [°C]', 'FontSize', 11);
    ylabel('h [W/(m^2·K)]', 'FontSize', 11);
    legend('Location', 'best', 'FontSize', 10);
    hold off;

    % 6.2 Lagrange 
    figure('Name', 'Interpolacja Lagrange''a');
    hold on; grid on;
    plot(dT_tab, h_tab, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Dane pomiarowe');
    plot(dT_dense, fh_lagrange(dT_dense), 'r', 'LineWidth', 2, 'DisplayName', 'Lagrange');
    title('Interpolacja Lagrange''a', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('\DeltaT [°C]', 'FontSize', 11);
    ylabel('h [W/(m^2·K)]', 'FontSize', 11);
    legend('Location', 'best', 'FontSize', 10);
    hold off;
    
    % 6.3 Splajn 
    figure('Name', 'Splajn kubiczny');
    hold on; grid on;
    plot(dT_tab, h_tab, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Dane pomiarowe');
    plot(dT_dense, fh_spline(dT_dense), 'm', 'LineWidth', 2, 'DisplayName', 'Splajn');
    title('Splajn kubiczny', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('\DeltaT [°C]', 'FontSize', 11);
    ylabel('h [W/(m^2·K)]', 'FontSize', 11);
    legend('Location', 'best', 'FontSize', 10);
    hold off;

    % 6.4 Wykres porównawczy wszystkich metod
    figure('Name', 'Porównanie metod');
    hold on; grid on;
    plot(dT_tab, h_tab, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Dane pomiarowe');
    plot(dT_dense, fh_mnk(dT_dense), 'Color', [1 0.5 0], 'LineWidth', 1.5, 'DisplayName', 'MNK');
    plot(dT_dense, fh_lagrange(dT_dense), 'r', 'LineWidth', 1.5, 'DisplayName', 'Lagrange');
    plot(dT_dense, fh_spline(dT_dense), 'm', 'LineWidth', 1.5, 'DisplayName', 'Splajn');
    title('Porównanie metod aproksymacji/interpolacji', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('\DeltaT [°C]', 'FontSize', 11);
    ylabel('h [W/(m^2·K)]', 'FontSize', 11);
    legend('Location', 'best', 'FontSize', 10);
    hold off;

    %% 7. Symulacja chłodzenia - parametry
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
    h_step = 0.001;

    %% 8. Symulacje dla każdej metody
    fprintf('\nRozpoczynanie symulacji...\n');
    
    [tM, xM] = metodaEuleraUlepszona_nielin(@f_nielin, t0, tK, h_step, x0, parametry, fh_mnk);
    fprintf('  - MNK: zakończono\n');
    
    [tL, xL] = metodaEuleraUlepszona_nielin(@f_nielin, t0, tK, h_step, x0, parametry, fh_lagrange);
    fprintf('  - Lagrange: zakończono\n');
    
    [tS, xS] = metodaEuleraUlepszona_nielin(@f_nielin, t0, tK, h_step, x0, parametry, fh_spline);
    fprintf('  - Splajn: zakończono\n\n');

    %% 9. Wykresy TB(t) i TW(t)
    
    % 9.1 Temperatura preta TB(t)
    figure('Name', 'Temperatura preta');
    hold on; grid on;
    plot(tM, xM(:,1), 'Color', [1 0.5 0], 'LineWidth', 1.8, 'DisplayName', 'MNK');
    plot(tL, xL(:,1), 'r', 'LineWidth', 1.8, 'DisplayName', 'Lagrange');
    plot(tS, xS(:,1), 'm', 'LineWidth', 1.8, 'DisplayName', 'Splajn');
    xlabel('t [s]', 'FontSize', 11); 
    ylabel('T_b [°C]', 'FontSize', 11);
    title('Temperatura preta T_b(t) - porównanie metod', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    hold off;

    % 9.2 Temperatura cieczy chlodzacej TW(t)
    figure('Name', 'Temperatura cieczy chlodzacej');
    hold on; grid on;
    plot(tM, xM(:,2), 'Color', [1 0.5 0], 'LineWidth', 1.8, 'DisplayName', 'MNK');
    plot(tL, xL(:,2), 'r', 'LineWidth', 1.8, 'DisplayName', 'Lagrange');
    plot(tS, xS(:,2), 'm', 'LineWidth', 1.8, 'DisplayName', 'Splajn');
    xlabel('t [s]', 'FontSize', 11); 
    ylabel('T_w [°C]', 'FontSize', 11);
    title('Temperatura cieczy chlodzacej T_w(t) - porównanie metod', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    hold off;

    %% 10. Obliczenie błędów bezwzględnych i względnych
    h_mnk_tab      = fh_mnk(dT_tab);
    h_lagrange_tab = fh_lagrange(dT_tab);
    h_spline_tab   = fh_spline(dT_tab);
    
    metody = {'MNK'; 'Lagrange'; 'Splajn'};
    H = [h_mnk_tab(:), h_lagrange_tab(:), h_spline_tab(:)];
    
    blad_bezwzgledny_sredni = zeros(3, 1);
    blad_bezwzgledny_max    = zeros(3, 1);
    blad_wzgledny_sredni    = zeros(3, 1);
    blad_wzgledny_max       = zeros(3, 1);
    
    for i = 1:3
        err_abs = abs(H(:,i) - h_tab(:));
        err_rel = abs((H(:,i) - h_tab(:)) ./ h_tab(:)) * 100;
    
        blad_bezwzgledny_sredni(i) = mean(err_abs);
        blad_bezwzgledny_max(i)    = max(err_abs);
        blad_wzgledny_sredni(i)    = mean(err_rel);
        blad_wzgledny_max(i)       = max(err_rel);
    end
    
    T_bledy = table( ...
        metody, ...
        blad_bezwzgledny_sredni, ...
        blad_bezwzgledny_max, ...
        blad_wzgledny_sredni, ...
        blad_wzgledny_max, ...
        'VariableNames', { ...
            'Metoda', ...
            'Blad_bezwzgl_sredni', ...
            'Blad_bezwzgl_max', ...
            'Blad_wzgl_sredni_proc', ...
            'Blad_wzgl_max_proc' ...
        } ...
    );

    fprintf('BŁĘDY BEZWZGLĘDNE I WZGLĘDNE - PODSUMOWANIE\n');
    disp(T_bledy);


end