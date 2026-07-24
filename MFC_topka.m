function [ek,e]=MFC_topka(at,bt,t,D, hp, d33,d31,dyhm,jdzb,jdzb1,dybh)
%----------------------------------------------------------------------------------------------   
% dyhm: node connectivity for each element  
% jdzb: coordinate of the nodes in the top surface of element
% jdzb1: coordinate of the nodes in the bottom surface of element
% dybh: the number of element whose stiffness matrix is calculated in this time of invoking. 
temp=0.577350292; % value needed by Gauss integral
gaosi=[ 1  1  1;
        1 -1  1;
       -1  1  1;
       -1 -1  1;
        1  1 -1;
        1 -1 -1;
       -1  1 -1;
       -1 -1 -1];

k=1.2;   

d2=[ d33 0 0 ;
    d31 0 0 ;
    0 0 0;
    0 0 0;
    0 0 0];
e=d2'*D;

% d2=[ 0 0 4e-10;
%     0 0 -1.7e-10;
%     0 0 0;
%     0 0 0;
%     0 0 0];
% e=0.36*d2'*D;

% d1=[ 0 0 4.67e-10; 
%     0 0 -2.1e-10;
%     0 0 0  ;
%     0 0 0;
%     0 0 0];
% e11=d1'*D;
% 
% d2=[ 0 0 -1.7e-10;
%     0 0 -1e-10;
%     0 0 0 ;
%     0 0 0;
%     0 0 0];
% e22=0.6*d2'*D;

Vp=[1/0.0005;
    0;
    0];

ek=zeros(40,1);
    for i=1:8
        for j=1:3
            v3i(i,j)=jdzb(dybh(dyhm,i),j)-jdzb1(dybh(dyhm,i),j);
        end
    end
    tempi=[1 0 0];
    for i=1:8
       temp1=cross(v3i(i,:)',tempi);
       xv2i(i,:)=temp1/norm(temp1);
       temp1=cross(xv2i(i,:),v3i(i,:));
       xv1i(i,:)=temp1/norm(temp1);
    end
    xv3i=v3i/t;
    
    for gaosii=1:8  % begin of Gauss integral
    zb1=temp*gaosi(gaosii,1);
    zb2=temp*gaosi(gaosii,2);
    zbs=temp*gaosi(gaosii,3);
    
    bb=bt;aa=at;
    zb3=((bb+aa+(bb-aa)*zbs)/2);

    ni(1)=1/4*(1+zb1)*(1+zb2)*(zb1+zb2-1);
    ni(2)=1/4*(1-zb1)*(1+zb2)*(-zb1+zb2-1);
    ni(3)=1/4*(1-zb1)*(1-zb2)*(-zb1-zb2-1);
    ni(4)=1/4*(1+zb1)*(1-zb2)*(zb1-zb2-1);
    ni(5)=1/2*(1-zb1^2)*(1+zb2);
    ni(6)=1/2*(1-zb1)*(1-zb2^2);
    ni(7)=1/2*(1-zb1^2)*(1-zb2);
    ni(8)=1/2*(1+zb1)*(1-zb2^2);
    v3=ni(1)*v3i(1,:)'+ni(2)*v3i(2,:)'+ni(3)*v3i(3,:)'+ni(4)*v3i(4,:)'+ni(5)*v3i(5,:)'+ni(6)*v3i(6,:)'+ni(7)*v3i(7,:)'+ni(8)*v3i(8,:)';
    temp1=cross(v3,tempi');
    tn2=sqrt(temp1(1)^2+temp1(2)^2+temp1(3)^2);
    xv2=temp1/tn2;
    temp1=cross(xv2,v3);
    tn1=sqrt(temp1(1)^2+temp1(2)^2+temp1(3)^2);
    xv1=temp1/tn1;
    
    tn3=sqrt(v3(1)^2+v3(2)^2+v3(3)^2);
    xv3=v3/tn3;
    theta=[xv1,xv2,xv3];
    
    for i=1:8
       for j=1:3
          zmtemp(i,j)=(jdzb1(dybh(dyhm,i),j)+jdzb(dybh(dyhm,i),j))/2;
       end
    end
    
    
    pni(1,1)=1/4*(1+zb2)*(zb1+zb2-1)+(1/4+1/4*zb1)*(1+zb2);
    pni(1,2)=(1/4+1/4*zb1)*(zb1+zb2-1)+(1/4+1/4*zb1)*(1+zb2);
    pni(1,3)=0;
    pni(2,1)=-1/4*(1+zb2)*(-zb1+zb2-1)-(1/4-1/4*zb1)*(1+zb2);
    pni(2,2)=(1/4-1/4*zb1)*(-zb1+zb2-1)+(1/4-1/4*zb1)*(1+zb2);
    pni(2,3)=0;
    pni(3,1)=-1/4*(1-zb2)*(-zb1-zb2-1)-(1/4-1/4*zb1)*(1-zb2);
    pni(3,2)=-(1/4-1/4*zb1)*(-zb1-zb2-1)-(1/4-1/4*zb1)*(1-zb2);
    pni(3,3)=0;
    pni(4,1)=1/4*(1-zb2)*(zb1-zb2-1)+(1/4+1/4*zb1)*(1-zb2);
    pni(4,2)=-(1/4+1/4*zb1)*(zb1-zb2-1)-(1/4+1/4*zb1)*(1-zb2);
    pni(4,3)=0;
    pni(5,1)=-zb1*(1+zb2);
    pni(5,2)=1/2-1/2*zb1^2;
    pni(5,3)=0;
    pni(6,1)=-1/2+1/2*zb2^2;
    pni(6,2)=-2*(1/2-1/2*zb1)*zb2;
    pni(6,3)=0;
    pni(7,1)=-zb1*(1-zb2);
    pni(7,2)=-1/2+1/2*zb1^2;
    pni(7,3)=0;
    pni(8,1)=1/2-1/2*zb2^2;
    pni(8,2)=-2*(1/2+1/2*zb1)*zb2;
    pni(8,3)=0;
    
    x1=zmtemp(1,1);
    y1=zmtemp(1,2);
    z1=zmtemp(1,3);
    x2=zmtemp(2,1);
    y2=zmtemp(2,2);
    z2=zmtemp(2,3);
    x3=zmtemp(3,1);
    y3=zmtemp(3,2);
    z3=zmtemp(3,3);
    x4=zmtemp(4,1);
    y4=zmtemp(4,2);
    z4=zmtemp(4,3);
    x5=zmtemp(5,1);
    y5=zmtemp(5,2);
    z5=zmtemp(5,3);
    x6=zmtemp(6,1);
    y6=zmtemp(6,2);
    z6=zmtemp(6,3);
    x7=zmtemp(7,1);
    y7=zmtemp(7,2);
    z7=zmtemp(7,3);
    x8=zmtemp(8,1);
    y8=zmtemp(8,2);
    z8=zmtemp(8,3);
    
    dx1=v3i(1,1);
    dy1=v3i(1,2);
    dz1=v3i(1,3);
    dx2=v3i(2,1);
    dy2=v3i(2,2);
    dz2=v3i(2,3);
    dx3=v3i(3,1);
    dy3=v3i(3,2);
    dz3=v3i(3,3);
    dx4=v3i(4,1);
    dy4=v3i(4,2);
    dz4=v3i(4,3);
    dx5=v3i(5,1);
    dy5=v3i(5,2);
    dz5=v3i(5,3);
    dx6=v3i(6,1);
    dy6=v3i(6,2);
    dz6=v3i(6,3);
    dx7=v3i(7,1);
    dy7=v3i(7,2);
    dz7=v3i(7,3);
    dx8=v3i(8,1);
    dy8=v3i(8,2);
    dz8=v3i(8,3);
     
    jtemp(1,1)=1/4*(1+zb2)*(zb1+zb2-1)*x1+(1/4+1/4*zb1)*(1+zb2)*x1-1/4*(1+zb2)*(-zb1+zb2-1)*x2-(1/4-1/4*zb1)*(1+zb2)*x2-1/4*(1-zb2)*(-zb1-zb2-1)*x3-(1/4-1/4*zb1)*(1-zb2)*x3+1/4*(1-zb2)*(zb1-zb2-1)*x4+(1/4+1/4*zb1)*(1-zb2)*x4-zb1*(1+zb2)*x5-1/2*(1-zb2^2)*x6-zb1*(1-zb2)*x7+1/2*(1-zb2^2)*x8+1/8*zb3*(1+zb2)*(zb1+zb2-1)*dx1+1/2*zb3*(1/4+1/4*zb1)*(1+zb2)*dx1-1/8*zb3*(1+zb2)*(-zb1+zb2-1)*dx2-1/2*zb3*(1/4-1/4*zb1)*(1+zb2)*dx2-1/8*zb3*(1-zb2)*(-zb1-zb2-1)*dx3-1/2*zb3*(1/4-1/4*zb1)*(1-zb2)*dx3+1/8*zb3*(1-zb2)*(zb1-zb2-1)*dx4+1/2*zb3*(1/4+1/4*zb1)*(1-zb2)*dx4-1/2*zb3*zb1*(1+zb2)*dx5-1/4*zb3*(1-zb2^2)*dx6-1/2*zb3*zb1*(1-zb2)*dx7+1/4*zb3*(1-zb2^2)*dx8;
    jtemp(1,2)=1/4*(1+zb2)*(zb1+zb2-1)*y1+(1/4+1/4*zb1)*(1+zb2)*y1-1/4*(1+zb2)*(-zb1+zb2-1)*y2-(1/4-1/4*zb1)*(1+zb2)*y2-1/4*(1-zb2)*(-zb1-zb2-1)*y3-(1/4-1/4*zb1)*(1-zb2)*y3+1/4*(1-zb2)*(zb1-zb2-1)*y4+(1/4+1/4*zb1)*(1-zb2)*y4-zb1*(1+zb2)*y5-1/2*(1-zb2^2)*y6-zb1*(1-zb2)*y7+1/2*(1-zb2^2)*y8+1/8*zb3*(1+zb2)*(zb1+zb2-1)*dy1+1/2*zb3*(1/4+1/4*zb1)*(1+zb2)*dy1-1/8*zb3*(1+zb2)*(-zb1+zb2-1)*dy2-1/2*zb3*(1/4-1/4*zb1)*(1+zb2)*dy2-1/8*zb3*(1-zb2)*(-zb1-zb2-1)*dy3-1/2*zb3*(1/4-1/4*zb1)*(1-zb2)*dy3+1/8*zb3*(1-zb2)*(zb1-zb2-1)*dy4+1/2*zb3*(1/4+1/4*zb1)*(1-zb2)*dy4-1/2*zb3*zb1*(1+zb2)*dy5-1/4*zb3*(1-zb2^2)*dy6-1/2*zb3*zb1*(1-zb2)*dy7+1/4*zb3*(1-zb2^2)*dy8;
    jtemp(1,3)=1/4*(1+zb2)*(zb1+zb2-1)*z1+(1/4+1/4*zb1)*(1+zb2)*z1-1/4*(1+zb2)*(-zb1+zb2-1)*z2-(1/4-1/4*zb1)*(1+zb2)*z2-1/4*(1-zb2)*(-zb1-zb2-1)*z3-(1/4-1/4*zb1)*(1-zb2)*z3+1/4*(1-zb2)*(zb1-zb2-1)*z4+(1/4+1/4*zb1)*(1-zb2)*z4-zb1*(1+zb2)*z5-1/2*(1-zb2^2)*z6-zb1*(1-zb2)*z7+1/2*(1-zb2^2)*z8+1/8*zb3*(1+zb2)*(zb1+zb2-1)*dz1+1/2*zb3*(1/4+1/4*zb1)*(1+zb2)*dz1-1/8*zb3*(1+zb2)*(-zb1+zb2-1)*dz2-1/2*zb3*(1/4-1/4*zb1)*(1+zb2)*dz2-1/8*zb3*(1-zb2)*(-zb1-zb2-1)*dz3-1/2*zb3*(1/4-1/4*zb1)*(1-zb2)*dz3+1/8*zb3*(1-zb2)*(zb1-zb2-1)*dz4+1/2*zb3*(1/4+1/4*zb1)*(1-zb2)*dz4-1/2*zb3*zb1*(1+zb2)*dz5-1/4*zb3*(1-zb2^2)*dz6-1/2*zb3*zb1*(1-zb2)*dz7+1/4*zb3*(1-zb2^2)*dz8;
    jtemp(2,1)=(1/4+1/4*zb1)*(zb1+zb2-1)*x1+(1/4+1/4*zb1)*(1+zb2)*x1+(1/4-1/4*zb1)*(-zb1+zb2-1)*x2+(1/4-1/4*zb1)*(1+zb2)*x2-(1/4-1/4*zb1)*(-zb1-zb2-1)*x3-(1/4-1/4*zb1)*(1-zb2)*x3-(1/4+1/4*zb1)*(zb1-zb2-1)*x4-(1/4+1/4*zb1)*(1-zb2)*x4+(1/2-1/2*zb1^2)*x5-2*(1/2-1/2*zb1)*zb2*x6-(1/2-1/2*zb1^2)*x7-2*(1/2+1/2*zb1)*zb2*x8+1/2*zb3*(1/4+1/4*zb1)*(zb1+zb2-1)*dx1+1/2*zb3*(1/4+1/4*zb1)*(1+zb2)*dx1+1/2*zb3*(1/4-1/4*zb1)*(-zb1+zb2-1)*dx2+1/2*zb3*(1/4-1/4*zb1)*(1+zb2)*dx2-1/2*zb3*(1/4-1/4*zb1)*(-zb1-zb2-1)*dx3-1/2*zb3*(1/4-1/4*zb1)*(1-zb2)*dx3-1/2*zb3*(1/4+1/4*zb1)*(zb1-zb2-1)*dx4-1/2*zb3*(1/4+1/4*zb1)*(1-zb2)*dx4+1/2*zb3*(1/2-1/2*zb1^2)*dx5-zb3*(1/2-1/2*zb1)*zb2*dx6-1/2*zb3*(1/2-1/2*zb1^2)*dx7-zb3*(1/2+1/2*zb1)*zb2*dx8;
    jtemp(2,2)=(1/4+1/4*zb1)*(zb1+zb2-1)*y1+(1/4+1/4*zb1)*(1+zb2)*y1+(1/4-1/4*zb1)*(-zb1+zb2-1)*y2+(1/4-1/4*zb1)*(1+zb2)*y2-(1/4-1/4*zb1)*(-zb1-zb2-1)*y3-(1/4-1/4*zb1)*(1-zb2)*y3-(1/4+1/4*zb1)*(zb1-zb2-1)*y4-(1/4+1/4*zb1)*(1-zb2)*y4+(1/2-1/2*zb1^2)*y5-2*(1/2-1/2*zb1)*zb2*y6-(1/2-1/2*zb1^2)*y7-2*(1/2+1/2*zb1)*zb2*y8+1/2*zb3*(1/4+1/4*zb1)*(zb1+zb2-1)*dy1+1/2*zb3*(1/4+1/4*zb1)*(1+zb2)*dy1+1/2*zb3*(1/4-1/4*zb1)*(-zb1+zb2-1)*dy2+1/2*zb3*(1/4-1/4*zb1)*(1+zb2)*dy2-1/2*zb3*(1/4-1/4*zb1)*(-zb1-zb2-1)*dy3-1/2*zb3*(1/4-1/4*zb1)*(1-zb2)*dy3-1/2*zb3*(1/4+1/4*zb1)*(zb1-zb2-1)*dy4-1/2*zb3*(1/4+1/4*zb1)*(1-zb2)*dy4+1/2*zb3*(1/2-1/2*zb1^2)*dy5-zb3*(1/2-1/2*zb1)*zb2*dy6-1/2*zb3*(1/2-1/2*zb1^2)*dy7-zb3*(1/2+1/2*zb1)*zb2*dy8;
    jtemp(2,3)=(1/4+1/4*zb1)*(zb1+zb2-1)*z1+(1/4+1/4*zb1)*(1+zb2)*z1+(1/4-1/4*zb1)*(-zb1+zb2-1)*z2+(1/4-1/4*zb1)*(1+zb2)*z2-(1/4-1/4*zb1)*(-zb1-zb2-1)*z3-(1/4-1/4*zb1)*(1-zb2)*z3-(1/4+1/4*zb1)*(zb1-zb2-1)*z4-(1/4+1/4*zb1)*(1-zb2)*z4+(1/2-1/2*zb1^2)*z5-2*(1/2-1/2*zb1)*zb2*z6-(1/2-1/2*zb1^2)*z7-2*(1/2+1/2*zb1)*zb2*z8+1/2*zb3*(1/4+1/4*zb1)*(zb1+zb2-1)*dz1+1/2*zb3*(1/4+1/4*zb1)*(1+zb2)*dz1+1/2*zb3*(1/4-1/4*zb1)*(-zb1+zb2-1)*dz2+1/2*zb3*(1/4-1/4*zb1)*(1+zb2)*dz2-1/2*zb3*(1/4-1/4*zb1)*(-zb1-zb2-1)*dz3-1/2*zb3*(1/4-1/4*zb1)*(1-zb2)*dz3-1/2*zb3*(1/4+1/4*zb1)*(zb1-zb2-1)*dz4-1/2*zb3*(1/4+1/4*zb1)*(1-zb2)*dz4+1/2*zb3*(1/2-1/2*zb1^2)*dz5-zb3*(1/2-1/2*zb1)*zb2*dz6-1/2*zb3*(1/2-1/2*zb1^2)*dz7-zb3*(1/2+1/2*zb1)*zb2*dz8;
    jtemp(3,1)=1/2*(1/4+1/4*zb1)*(1+zb2)*(zb1+zb2-1)*dx1+1/2*(1/4-1/4*zb1)*(1+zb2)*(-zb1+zb2-1)*dx2+1/2*(1/4-1/4*zb1)*(1-zb2)*(-zb1-zb2-1)*dx3+1/2*(1/4+1/4*zb1)*(1-zb2)*(zb1-zb2-1)*dx4+1/2*(1/2-1/2*zb1^2)*(1+zb2)*dx5+1/2*(1/2-1/2*zb1)*(1-zb2^2)*dx6+1/2*(1/2-1/2*zb1^2)*(1-zb2)*dx7+1/2*(1/2+1/2*zb1)*(1-zb2^2)*dx8;
    jtemp(3,2)=1/2*(1/4+1/4*zb1)*(1+zb2)*(zb1+zb2-1)*dy1+1/2*(1/4-1/4*zb1)*(1+zb2)*(-zb1+zb2-1)*dy2+1/2*(1/4-1/4*zb1)*(1-zb2)*(-zb1-zb2-1)*dy3+1/2*(1/4+1/4*zb1)*(1-zb2)*(zb1-zb2-1)*dy4+1/2*(1/2-1/2*zb1^2)*(1+zb2)*dy5+1/2*(1/2-1/2*zb1)*(1-zb2^2)*dy6+1/2*(1/2-1/2*zb1^2)*(1-zb2)*dy7+1/2*(1/2+1/2*zb1)*(1-zb2^2)*dy8;
    jtemp(3,3)=1/2*(1/4+1/4*zb1)*(1+zb2)*(zb1+zb2-1)*dz1+1/2*(1/4-1/4*zb1)*(1+zb2)*(-zb1+zb2-1)*dz2+1/2*(1/4-1/4*zb1)*(1-zb2)*(-zb1-zb2-1)*dz3+1/2*(1/4+1/4*zb1)*(1-zb2)*(zb1-zb2-1)*dz4+1/2*(1/2-1/2*zb1^2)*(1+zb2)*dz5+1/2*(1/2-1/2*zb1)*(1-zb2^2)*dz6+1/2*(1/2-1/2*zb1^2)*(1-zb2)*dz7+1/2*(1/2+1/2*zb1)*(1-zb2^2)*dz8;
    invj=inv(jtemp); 
    dni=invj*[pni(1,1) pni(2,1) pni(3,1) pni(4,1) pni(5,1) pni(6,1) pni(7,1) pni(8,1);
                    pni(1,2) pni(2,2) pni(3,2) pni(4,2) pni(5,2) pni(6,2) pni(7,2) pni(8,2);
                    0 0 0 0 0 0 0 0];
    dmi=invj*[zb3*pni(1,1) zb3*pni(2,1) zb3*pni(3,1) zb3*pni(4,1) zb3*pni(5,1) zb3*pni(6,1) zb3*pni(7,1) zb3*pni(8,1);
              zb3*pni(1,2) zb3*pni(2,2) zb3*pni(3,2) zb3*pni(4,2) zb3*pni(5,2) zb3*pni(6,2) zb3*pni(7,2) zb3*pni(8,2);
              ni(1) ni(2) ni(3) ni(4) ni(5) ni(6) ni(7) ni(8)];
          
    for bi=1:8
        for s=1:3
        alpha(s)=theta(1,s)*dni(1,bi)+theta(2,s)*dni(2,bi)+theta(3,s)*dni(3,bi);
        beta(s)=(theta(1,s)*dmi(1,bi)+theta(2,s)*dmi(2,bi)+theta(3,s)*dmi(3,bi))*t/2;
        gama(s)=theta(1,s)*xv1i(bi,1)+theta(2,s)*xv1i(bi,2)+theta(3,s)*xv1i(bi,3);
        lmd(s)=theta(1,s)*xv2i(bi,1)+theta(2,s)*xv2i(bi,2)+theta(3,s)*xv2i(bi,3);
        end
        btemp(:,:,bi)=[theta(1,1)*alpha(1) theta(2,1)*alpha(1) theta(3,1)*alpha(1) beta(1)*gama(1) beta(1)*lmd(1);
                       theta(1,2)*alpha(2) theta(2,2)*alpha(2) theta(3,2)*alpha(2) beta(2)*gama(2) beta(2)*lmd(2);
                       theta(1,1)*alpha(2)+theta(1,2)*alpha(1) theta(2,1)*alpha(2)+theta(2,2)*alpha(1) theta(3,1)*alpha(2)+theta(3,2)*alpha(1) beta(1)*gama(2)+beta(2)*gama(1) beta(1)*lmd(2)+beta(2)*lmd(1);
                       theta(1,2)*alpha(3)+theta(1,3)*alpha(2) theta(2,2)*alpha(3)+theta(2,3)*alpha(2) theta(3,2)*alpha(3)+theta(3,3)*alpha(2) beta(2)*gama(3)+beta(3)*gama(2) beta(2)*lmd(3)+beta(3)*lmd(2);
                       theta(1,3)*alpha(1)+theta(1,1)*alpha(3) theta(2,3)*alpha(1)+theta(2,1)*alpha(3) theta(3,3)*alpha(1)+theta(3,1)*alpha(3) beta(3)*gama(1)+beta(1)*gama(3) beta(3)*lmd(1)+beta(1)*lmd(3)];
    end
    b=[btemp(:,:,1) btemp(:,:,2) btemp(:,:,3) btemp(:,:,4) btemp(:,:,5) btemp(:,:,6) btemp(:,:,7) btemp(:,:,8)];
    djt=det(jtemp);
    tk=b'*e'*Vp*djt;
    ek=ek+0.5*(bb-aa)*tk;
    end