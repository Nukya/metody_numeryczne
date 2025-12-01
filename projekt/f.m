%% Równania stanu


function dx = f(~, x, parametry)
    Tb = x(1);
    Tw = x(2);
    
    wspWnikania = parametry.h;
    powierzchnia = parametry.A;
    masaPreta = parametry.mb;
    masaOleju = parametry.mw;
    cieploPreta = parametry.cb;
    cieploOleju = parametry.cw;
    
    dTb = -(wspWnikania*powierzchnia)/(masaPreta*cieploPreta) * (Tb - Tw);
    dTw =  +(wspWnikania*powierzchnia)/(masaOleju*cieploOleju) * (Tb - Tw);
    
    dx = [dTb; dTw];
end