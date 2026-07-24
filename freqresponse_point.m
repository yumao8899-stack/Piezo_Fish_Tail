function Hd = freqresponse_point(K, M, C, F, Omega, nummode, ith)
%-------------------------------------------------------------------------
% Purpose:
% Calculate frequency response at specific response point for structural system
% with force vector input
% 
% Synopsis:
%       Hd = freqresponse_point(K, M, C, F, Omega, nummode, ith)
%
% Variable Description:
%       Input parameters
%       K - System stiffness matrix (nsdof x nsdof)
%       M - System mass matrix (nsdof x nsdof)  
%       C - System damping matrix (nsdof x nsdof)
%       F - Force vector (nsdof x 1) or (nsdof x n) for frequency-dependent forces
%       Omega - Frequency range vector (1 x n)
%       nummode - number of extracted modes
%       ith - response degree of freedom
%
%       Output parameters
%       Hd - Displacement frequency response at ith DOF (n x 1)
%
% Author: Adapted from Dr. XU Bin's code
%--------------------------------------------------------------------------

% Input validation
if nargin ~= 7
    error('Incorrect number of input arguments. Expected 7 inputs.')
end

nsdof = length(K);
n = length(Omega);

% Check force vector dimensions
if size(F, 1) ~= nsdof
    error('Force vector F must have the same number of rows as system DOF');
end

if size(F, 2) == 1
    % Constant force vector across all frequencies
    F = repmat(F, 1, n);
elseif size(F, 2) ~= n
    error('Force vector F must have either 1 column or same number of columns as Omega');
end

% Solve eigenvalue problem
nm=20;
[V,d]=eigs(k,m,nm,'SM');
tempd = diag(d);
[~, sortindex] = sort(tempd);
V = V(:, sortindex);

% Normalize modes
Factor = diag(V' * M * V);
Vnorm = V * diag(1./sqrt(Factor));

% Project damping matrix to modal coordinates
Dampr = Vnorm' * C * Vnorm;

% Initialize output
Hd = zeros(n, 1);

% Calculate frequency response using modal superposition
for p = 1:n
    Hd_temp = 0;
    for q = 1:nummode
        % Modal contribution to response at ith DOF
        modal_response = Vnorm(ith, q) * (Vnorm(:, q)' * F(:, p)) / ...
                        (lambda(q) - Omega(p)^2 + 1i * Dampr(q, q) * Omega(p));
        Hd_temp = Hd_temp + modal_response;
    end
    Hd(p) = Hd_temp;
end

end