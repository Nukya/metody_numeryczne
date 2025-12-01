function projekt_chlodzenie_cz4()
    clc; close all;

    %% === 1. Model h(ΔT) — splajn ===
    [dT_tab, h_tab] = dane_h();
    N = numel(dT_tab);
    dT_eq = linspace(min(dT_tab), max(dT_tab), N);
    h_eq  = interpolacja_liniowa(dT_tab, h_tab, dT_eq);

    alpha = 0; beta = 0;
    S = funkcje_sklejane(dT_eq, h_eq, alpha, beta);
    fh = @(dT) S.eval(dT);     % model h(ΔT)

    %% === 2. Parametry fizyczne ===
    parametry.A  = 0.0109;
    parametry.mb = 0.25;
    parametry.cb = 0.29;
    parametry.cw = 4.1813;

    Tb0 = 1200;
    Tw0 = 25;
    tK = 0.7;
    dt = 0.001;

    %% === 3. Funkcja Tb(mw) — UŻYWA nieliniowego modelu ===
    Tb_po07 = @(mw) symuluj_nielin(Tb0, Tw0, mw, tK, parametry, dt, fh);

    %% === 4. Równanie F(mw) = Tb(mw) - 125 ===
    T_cel = 125;
    Tb_po07 = @(mw) symuluj_nielin(Tb0, Tw0, mw, tK, parametry, dt, fh);

    %% === 5. Parametry Newtona ===
    mw0 = 2.0;
    d_mw = 0.02;
    tol = 1e-3;
    maxIter = 20;

    %% === 6. Iteracje Newtona ===
    mw_opt = newton_raphson(Tb_po07, mw0, d_mw, tol, maxIter, T_cel);

    fprintf("\n>>> Minimalna masa oleju: mw ≈ %.4f kg\n", mw_opt);
end
