function Tb_end = symuluj_nielin(Tb0, Tw0, mw, tK, parametry, h, fh)
    p = parametry;
    p.mw = mw;

    x0 = [Tb0 ; Tw0];

    [t, x] = metodaEuleraUlepszona_nielin(@f_nielin, 0, tK, h, x0, p, fh);
    Tb_end = x(end,1);
end
