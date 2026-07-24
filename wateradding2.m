function [Mw_modal, Cw_modal] =  wateradding2(Phi, zmtemp, rho_fluid, mu, omega, a1, a2, b1)

nNode   = 8;
nDOF    = size(Phi,1);
nModal  = size(Phi,2);

Mw_modal = zeros(nModal, nModal);
Cw_modal = zeros(nModal, nModal);

% ---------------- Gaussian quadrature ----------------
gp = [-0.774596669241483 0 0.774596669241483];
gw = [0.555555555555556 0.888888888888889 0.555555555555556];

% =========================================================
% LOOP over Gauss points
% =========================================================
for i = 1:3
    for j = 1:3

        xi = gp(i);
        eta = gp(j);
        W = gw(i)*gw(j);

        % ---------------- shape functions ----------------
        [N, dN] = shape8(xi, eta);

        % ---------------- Jacobian + geometry ----------------
        J11=0; J12=0; J21=0; J22=0;
        Xg=0;

        for k=1:nNode
            J11 = J11 + dN(k,1)*zmtemp(k,1);
            J12 = J12 + dN(k,1)*zmtemp(k,2);
            J21 = J21 + dN(k,2)*zmtemp(k,1);
            J22 = J22 + dN(k,2)*zmtemp(k,2);
            Xg  = Xg  + N(k)*zmtemp(k,1);
        end

        dA = J11*J22 - J12*J21;

        % ---------------- section properties ----------------
        b_local = get_fin_width(Xg);
        b_gap   = get_gap_width(Xg);

        beta  = (rho_fluid * omega * b_local^2) / (2*pi*mu);
        delta = b_gap / b_local;

        ReTheta = 1.02 + a1*delta + (2.45 + a2*delta)*beta^(-0.5);
        ImTheta = 2.49 + b1*delta*beta^(-0.5);

        m_area = (pi/4) * rho_fluid * b_local * ReTheta;
        c_area = (pi/4) * rho_fluid * b_local * ImTheta * omega;

        % =====================================================
        %  模态截断
        % =====================================================
        for a = 1:nModal
            for b = 1:nModal

                Phi_a = 0;
                Phi_b = 0;

                for k = 1:nNode
                    idx = (k-1)*5 + 3; % w DOF
                    Phi_a = Phi_a + Phi(idx,a) * N(k);
                    Phi_b = Phi_b + Phi(idx,b) * N(k);
                end

                Mw_modal(a,b) = Mw_modal(a,b) + m_area * Phi_a * Phi_b * dA * W;
                Cw_modal(a,b) = Cw_modal(a,b) + c_area * Phi_a * Phi_b * dA * W;

            end
        end

    end
end
end

function [N, dN] = shape8(xi, eta)
N(1)=1/4*(1+xi)*(1+eta)*(xi+eta-1);
N(2)=1/4*(1-xi)*(1+eta)*(-xi+eta-1);
N(3)=1/4*(1-xi)*(1-eta)*(-xi-eta-1);
N(4)=1/4*(1+xi)*(1-eta)*(xi-eta-1);
N(5)=1/2*(1-xi^2)*(1+eta);
N(6)=1/2*(1-xi)*(1-eta^2);
N(7)=1/2*(1-xi^2)*(1-eta);
N(8)=1/2*(1+xi)*(1-eta^2);
dN=zeros(8,2);

% d/dxi
dN(1,1)=1/4*(1+eta)*(2*xi+eta);
dN(2,1)=1/4*(1+eta)*(-2*xi+eta);
dN(3,1)=1/4*(1-eta)*(-2*xi-eta);
dN(4,1)=1/4*(1-eta)*(2*xi-eta);
dN(5,1)=-xi*(1+eta);
dN(6,1)=-1/2*(1-eta^2);
dN(7,1)=-xi*(1-eta);
dN(8,1)= 1/2*(1-eta^2);

% d/deta
dN(1,2)=1/4*(1+xi)*(xi+2*eta);
dN(2,2)=1/4*(1-xi)*( -xi+2*eta);
dN(3,2)=1/4*(1-xi)*( -xi-2*eta);
dN(4,2)=1/4*(1+xi)*( xi-2*eta);
dN(5,2)=1/2*(1-xi^2);
dN(6,2)=-eta*(1-xi);
dN(7,2)=-1/2*(1-xi^2);
dN(8,2)=-eta*(1+xi);
end

function b_local = get_fin_width(x)
bm = 0.049; be = 0.040;
lm = 0.020; le = 0.055;
if x <= bm
    b_local = lm;
else
    t = min(max((x - bm)/be, 0), 1);
    S = 3*t^2 - 2*t^3;
    b_local = lm + (le - lm)*S;
end
end



function b_gap = get_gap_width(x)
L = 0.089;
df = 0.010;
le = 0.055;
x0 = L - df;
if x <= x0
    b_gap = 0;
else
    x = min(x, L);
    b_gap = le * sqrt(max(0, 1 - (L - x)/df));
end
end