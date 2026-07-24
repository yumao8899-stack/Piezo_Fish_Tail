function Q = Q_gexiangyixing(E1, E2, G23, G12, G13, nu12)
    k = 6/5;
    nu21 = (E2/E1)*nu12;
    den = 1 - nu12*nu21;
    Q = zeros(5);
    Q(1,1) = E1 / den;
    Q(1,2) = nu12 * E2 / den;
    Q(2,1) = Q(1,2);
    Q(2,2) = E2 / den;
    Q(3,3) = G23;
    Q(4,4) = G12 / k;
    Q(5,5) = G13 / k;
end