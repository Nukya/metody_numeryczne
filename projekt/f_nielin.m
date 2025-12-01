function dx = f_nielin(~, x, parametry, fh)
    Tb = x(1);
    Tw = x(2);

    dT = Tb - Tw;
    h_dyn = fh(dT);

    A  = parametry.A;
    mb = parametry.mb;
    mw = parametry.mw;
    cb = parametry.cb;
    cw = parametry.cw;

    dTb = -(h_dyn*A)/(mb*cb) * (Tb - Tw);
    dTw =  +(h_dyn*A)/(mw*cw) * (Tb - Tw);

    dx = [dTb; dTw];
end
