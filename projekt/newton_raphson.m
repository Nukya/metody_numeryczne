function mw_opt = newton_raphson(Tb_fun, mw0, d_mw, tol, maxIter, Tcel)
%--------------------------------------------------------------------------
%
% mw_{i+1} = mw_i - ( Tb(mw_i) - T_cel ) / Tb'(mw_i)
%
% gdzie
%
% Tb'(mw) ≈ [ Tb(mw + Δmw) – Tb(mw) ] / Δmw
%
%--------------------------------------------------------------------------
    
    mw = mw0;

    fprintf("Iteracje Newtona:\n");

    for i = 1:maxIter
        
        Tb_val  = Tb_fun(mw);          % Tb(mw)
        Tb_plus = Tb_fun(mw + d_mw);   % Tb(mw + Δmw)

        % pochodna numeryczna
        Tb_prime = (Tb_plus - Tb_val) / d_mw;
        
        deltaTb = Tb_val - Tcel;

        fprintf(" i=%2d: mw = %.5f kg, Tb = %.4f, Tb' = %.6f, deltaTb = %.6f \n", ...
                 i, mw, Tb_val, Tb_prime, deltaTb);

        % warunek stopu
        if abs(deltaTb) < tol
            break;
        end

        if abs(Tb_prime) < 1e-10
            warning("Pochodna bliska zeru — iteracje przerwane.");
            break;
        end

        mw_new = mw - (Tb_val - Tcel) / Tb_prime;

        if mw_new <= 0
            warning("Masz masa ujemna — korekta.");
            mw_new = mw / 2;
        end

        mw = mw_new;
    end

    mw_opt = mw;
end
