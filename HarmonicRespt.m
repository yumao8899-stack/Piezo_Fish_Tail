function [eta,y]=HarmonicRespt(kk,mm,fd,omega0,t,C,q0,dq0,CC,omega,V)
t=t';
nn=20;
[nstep,n2]=size(t);
V1=[V(:,1:nn)];                                     % truncate the modal vectors
Factor=diag(V1'*mm*V1);
Vnorm=V1*inv(sqrt(diag(Factor)));                    %  eigenvectors are normalized
omega2=diag(sqrt(Vnorm'*kk*Vnorm));                         % natural frequencies
           
Fnorm=Vnorm'*fd;                                      % modal input force vector
Modamp=Vnorm'*(CC)*Vnorm;                            % form the Rayleigh damping
zeta=diag((1/2)*Modamp*inv(diag(omega2)));  

eta0=Vnorm'*mm*q0; deta0=Vnorm'*mm*dq0;
eta=zeros(nstep,nn);

phase0=omega0*t;

for i=1:nn                                  % responses are obtained for n modes
  gama=omega0/omega(i);
  omegad=omega(i)*sqrt(1-zeta(i)^2);
  phase=omegad*t;
  Exx=exp(-zeta(i)*omega(i)*t);
  C1=eta0(i);
  C2=(deta0(i)+eta0(i)*zeta(i)*omega(i))/omegad;
  X0=sqrt((1-gama^2)^2+(2*zeta(i)*gama)^2);
  XX=Fnorm(i)/(omega(i)^2*X0);%
  XP=atan((2*zeta(i)*gama)/(1-gama^2));
  D1=(zeta(i)*omega(i)*cos(XP)+omega0*sin(XP))/omegad;
  D2=cos(XP);
  eta(:,i)=C1*Exx.*cos(phase)+C2*Exx.*sin(phase)...
           -XX*Exx.*(D1*sin(phase)+D2*cos(phase))...
           +XX*cos(phase0-XP);
end
eta=eta';
y=C*Vnorm*eta;
 

