function S = funkcje_sklejane(xi, yi, alpha, beta)
% Wejście:
%   xi    – węzły (równomierne)
%   yi    – wartości funkcji
%   alpha – pochodna w pierwszym węźle
%   beta  – pochodna w ostatnim węźle
%
% Wyjście:
%   S – struktura
%       S.xi    – rozszerzone węzły
%       S.wspC  – współczynniki c_i
%       S.h     – krok h
%       S.eval  – uchwyt funkcji S.eval(x)

    % --- parametry ---
    a = xi(1);
    b = xi(end);
    n = length(xi) - 1;
    h = (b - a) / n;


    % 1. Budowa macierzy A (trójdiagonalnej)

    A = zeros(n+1);
    A(1,1) = 4;  
    A(1,2) = 2;

    for i = 2:n
        A(i,i-1) = 1;
        A(i,i)   = 4;
        A(i,i+1) = 1;
    end

    A(n+1, n) = 2; 
    A(n+1,n+1) = 4;


    % 2. Wektor prawej strony d

    d = zeros(n+1,1);
    d(1) = yi(1) + h/3 * alpha;
    d(2:n) = yi(2:n);
    d(n+1) = yi(end) - h/3 * beta;


    % 3. Rozwiązanie A * C = d

    C = A \ d;

    % Rozszerzenie współczynników c o c(-1) i c(n+1)
    c_m1  = C(2) - h/3 * alpha;
    c_np1 = C(end-1) + h/3 * beta;

    wspC = [c_m1; C; c_np1];

    % Rozszerzenie węzłów
    xi_ext = [xi(1)-h, xi(:)', xi(end)+h];


    % 4. Zapisanie do struktury

    S.xi = xi_ext;
    S.wspC = wspC;
    S.h = h;

    % uchwyt — woła lokalną funkcję S_eval()
    S.eval = @(x) S_eval(S, x);

    % KONIEC GŁÓWNEJ FUNKCJI
end


function val = S_eval(S, x)
% S_eval — liczy wartość splajnu S(x)

    xi = S.xi;
    C  = S.wspC;
    h  = S.h;

    val = 0;
    for k = 1:length(xi)
        val = val + C(k) * phi_local(xi(k), h, x);
    end
end


function result = phi_local(xi, h, x)
    % wektorowe wejście
    x = x(:)';  

    xi_m2 = xi - 2*h;
    xi_m1 = xi - h;
    xi_p1 = xi + h;
    xi_p2 = xi + 2*h;

    y = zeros(size(x));

    idx = (x >= xi_m2 & x <= xi_m1);
    y(idx) = (x(idx) - xi_m2).^3;

    idx = (x >= xi_m1 & x <= xi);
    y(idx) = (x(idx) - xi_m2).^3 - 4*(x(idx) - xi_m1).^3;

    idx = (x >= xi & x <= xi_p1);
    y(idx) = (xi_p2 - x(idx)).^3 - 4*(xi_p1 - x(idx)).^3;

    idx = (x >= xi_p1 & x <= xi_p2);
    y(idx) = (xi_p2 - x(idx)).^3;

    result = y / h^3;
end

