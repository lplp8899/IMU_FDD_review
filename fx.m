function F=fx(x,u,g)
    phi=x(4);theta=x(5);psi=x(6);
    axm=u(1);aym=u(1);azm=u(1);pm=u(1);qm=u(1);rm=u(1);
    F(1,1) = axm*cos(theta)*cos(psi) + aym*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi)) + azm*(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi));
    F(2,1) = axm*cos(theta)*sin(psi) + aym*(sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi)) + azm*(cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi));
    F(3,1) = g-axm*sin(theta) + aym*sin(phi)*cos(theta) + azm*cos(phi)*cos(theta);
    F(4,1) = pm + qm*sin(phi)*tan(theta) + rm*cos(phi)*tan(theta);
    F(5,1) = qm*cos(phi) - rm*sin(phi);
    F(6,1) = qm*sin(phi)*sec(theta) + rm*cos(phi)*sec(theta);
    
    
    
%     BodyToEarth = [ cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
%                     cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
%                    -sin(theta)          sin(phi)*cos(theta)                            cos(phi)*cos(theta)];
%                
%     Eular_Lin =  [  1  sin(phi)*tan(theta)  cos(phi)*tan(theta);
%                     0  cos(phi)            -sin(phi);
%                     0  sin(phi)*sec(theta)  cos(phi)*sec(theta)];
%                            
%     G1_con = [ BodyToEarth zeros(3,3);
%                zeros(3,3) Eular_Lin];

end