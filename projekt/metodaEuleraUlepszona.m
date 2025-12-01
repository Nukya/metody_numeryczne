%% Ulepszona metoda Eulera

function [t, x] = metodaEuleraUlepszona(fhandle, t0, tK, h, x0, parametry)
    N = floor((tK - t0)/h) + 1;
    t = linspace(t0, t0 + h*(N-1), N).';
    x = zeros(N, length(x0));
    x(1,:) = x0.';
    
    for k = 1:N-1
        fk1 = fhandle(t(k), x(k,:).', parametry);    
        x_temp = x(k,:).' + h * fk1;                     
        fk2 = fhandle(t(k)+h, x_temp, parametry);         
    
        x(k+1,:) = x(k,:) + (h/2) * (fk1.' + fk2.');   
    end
end