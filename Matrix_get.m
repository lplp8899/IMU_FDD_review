function [G1_con,G2_con] = Matrix_get(x)

    phi = x(4);
    theta = x(5);
    psi = x(6);
    
    BodyToEarth = [ cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
                    cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
                   -sin(theta)          sin(phi)*cos(theta)                            cos(phi)*cos(theta)];
               
    Eular_Lin =  [  1  sin(phi)*tan(theta)  cos(phi)*tan(theta);
                    0  cos(phi)            -sin(phi);
                    0  sin(phi)*sec(theta)  cos(phi)*sec(theta)];
                           
    G1_con = [ BodyToEarth zeros(3,3);
               zeros(3,3) Eular_Lin];
          
    G2_con = -G1_con;
                   
end
