function s = lagrange(xp, yp, x)
    n = length(xp);
    s = zeros(size(x));
    for k = 1:n
        Lk = ones(size(x));
        for i = 1:n
            if i ~= k
                Lk = Lk .* (x - xp(i)) ./ (xp(k) - xp(i));
            end
        end
        s = s + yp(k) * Lk;
    end
end
