function y = wartosc_wielomianu(wsp, x)
    wsp = wsp(:).';
    n = length(wsp);
    y = zeros(size(x));

    for k = 1:numel(x)
        s = 0;
        for i = 1:n
            s = s + wsp(i) * x(k)^(n-i);
        end
        y(k) = s;
    end
end
