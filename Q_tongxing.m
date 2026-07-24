function Q = Q_tongxing(E, nu)
    k = 6/5;
    G = E / (2*(1 + nu));
    Q = zeros(5);
    Q(1,1) = E / (1 - nu^2);
    Q(1,2) = nu * Q(1,1);
    Q(2,1) = Q(1,2);
    Q(2,2) = Q(1,1);
    Q(3,3) = G;
    Q(4,4) = G / k;
    Q(5,5) = G / k;
end