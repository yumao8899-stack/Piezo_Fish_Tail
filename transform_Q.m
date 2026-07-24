function Qr = transform_Q(Q, theta)
theta=deg2rad(theta);
m = cos(theta); n = sin(theta);
T1 = [m^2, n^2, 2*m*n;
    n^2, m^2, -2*m*n;
    -m*n, m*n, m^2 - n^2];
T2 = [m,-n;n,m];
T6 = blkdiag(T1, T2);
Qr = T6 * Q * T6';
end