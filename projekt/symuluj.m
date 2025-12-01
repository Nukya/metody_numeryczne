%% Symulacja do t_obs

function Tb_koniec = symuluj(Tb0,Tw0,mw,czasObserwacji,parametry,krokCzasu)
    parametryLokalne = parametry;
    parametryLokalne.mw = mw;
    
    [~, rozw] = metodaEuleraUlepszona(@f, 0, czasObserwacji, krokCzasu, [Tb0;Tw0], parametryLokalne);
    Tb_koniec = rozw(end,1);
end
