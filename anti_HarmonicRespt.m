function [eta,y]=anti_HarmonicRespt(kk,mm,fd,omega0,t,C,q0,dq0,CC,omega,V,phi)
t=t';
nn=5;
[nstep,~]=size(t);
[~, n_forces] = size(fd);  % 获取激励力的数量

V1=[V(:,1:nn)];                                     
Factor=diag(V1'*mm*V1);
Vnorm=V1*inv(sqrt(diag(Factor)));                    
omega2=diag(sqrt(Vnorm'*kk*Vnorm));                         
           
Modamp=Vnorm'*(CC)*Vnorm;                            
zeta=diag((1/2)*Modamp*inv(diag(omega2)));  

eta0=Vnorm'*mm*q0; 
deta0=Vnorm'*mm*dq0;

eta = zeros(nstep, nn);  

% 先计算初始条件引起的自由振动响应
for i=1:nn
  omegad = omega(i)*sqrt(1-zeta(i)^2);
  phase = omegad*t ;
  Exx = exp(-zeta(i)*omega(i)*t);
  C1 = eta0(i);
  C2 = (deta0(i)+eta0(i)*zeta(i)*omega(i))/omegad;
  eta(:,i) = (C1*Exx.*cos(phase) + C2*Exx.*sin(phase))';  
end

% 对每个激励力分别计算强迫响应，然后叠加
for j=1:n_forces                                    
  Fnorm = Vnorm'*fd(:,j);                          
  phase0 = omega0(j)*t+ phi(j,1);                            
  
  for i=1:nn                                        
    gama = omega0(j)/omega(i);
    omegad = omega(i)*sqrt(1-zeta(i)^2);
    phase = omegad*t ;
    Exx = exp(-zeta(i)*omega(i)*t);
    
    X0 = sqrt((1-gama^2)^2+(2*zeta(i)*gama)^2);
    XX = Fnorm(i)/(omega(i)^2*X0);
    XP = atan((2*zeta(i)*gama)/(1-gama^2));
    D1 = (zeta(i)*omega(i)*cos(XP)+omega0(j)*sin(XP))/omegad;
    D2 = cos(XP);
    
    % 累加强迫响应
    forced_response = -XX*Exx.*(D1*sin(phase)+D2*cos(phase)) ...
                      + XX*cos(phase0-XP);
    
    eta(:,i) = eta(:,i) + forced_response;
  end
end

eta = eta';  
y = C*Vnorm*eta;