function projekt_chlodzenie_cz4()
    clc; close all;
    fprintf('CZĘŚĆ 4: OPTYMALIZACJA MASY OLEJU\n\n');
    
    %%  1. Model h(ΔT) — MNK (najlepsza z Części 3)
    [dT_tab, h_tab] = dane_h();
    stopien = 5;
    [wsp_mnk, ~] = aproksymacja_wielomianowa(dT_tab, h_tab, stopien);
    fh = @(dT) wartosc_wielomianu(wsp_mnk, dT);

    %%  2. Parametry
    parametry.A  = 0.0109;
    parametry.mb = 0.25;
    parametry.cb = 0.29;
    parametry.cw = 4.1813;
    Tb0 = 1200;
    Tw0 = 25;
    tK = 0.7;
    dt = 0.001;
    T_cel = 125;
    %%  3. Funkcja Tb(mw)
    Tb_po07 = @(mw) symuluj_nielin(Tb0, Tw0, mw, tK, parametry, dt, fh);

    %%  4. Parametry Newtona
    mw0 = 2.0;
    d_mw = 0.02;
    tol = 1e-3;
    maxIter = 20;

    %%  5. Metoda Newtona
    mw_opt = newton_raphson(Tb_po07, mw0, d_mw, tol, maxIter, T_cel);

    %%  6. Weryfikacja
    Tb_final = Tb_po07(mw_opt);
    blad = abs(Tb_final - T_cel);
    fprintf('\n========================================\n');
    fprintf('WYNIKI:\n\n');
    fprintf('Minimalna masa oleju: mw = %.4f kg\n', mw_opt);
    fprintf('Temperatura końcowa: %.6f°C\n', Tb_final);
    fprintf('Temperatura docelowa: %.0f°C\n', T_cel);
    fprintf('Błąd: %.6f°C\n', blad);

    %%  7. Wykres rozwiązania
    parametry.mw = mw_opt;
    x0 = [Tb0; Tw0];
    [t, x] = metodaEuleraUlepszona_nielin(@f_nielin, 0, tK, dt, x0, parametry, fh);
    figure('Name', 'optymalne rozwiazanie temperatury');
    subplot(2,1,1);
        plot(t, x(:,1), 'b', 'LineWidth', 2);
        hold on;
        plot(tK, T_cel, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        xlabel('t [s]');
        ylabel('T_b [°C]');
        title(sprintf('Temperatura pręta (m_w = %.4f kg)', mw_opt));
        legend('T_b(t)', sprintf('Cel: %.0f°C', T_cel));
        grid on;
    subplot(2,1,2);
        plot(t, x(:,2), 'r', 'LineWidth', 2);
        xlabel('t [s]');
        ylabel('T_w [°C]');
        title('Temperatura oleju');
        grid on;
end