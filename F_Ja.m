function F=F_Ja(x,uk,gk)

[num,~]  =  size(x);

UN = x(4-3,:); UE = x(5-3,:); UD = x(6-3,:);
phi = x(7-3,:); theta = x(8-3,:); psi = x(9-3,:);

if num  == 15 || num  == 12
    la_x = x(10-3,:); la_y = x(11-3,:); la_z = x(12-3,:);
    la_p = x(13-3,:); la_q = x(14-3,:); la_r = x(15-3,:);
%     uwind = x(13,:);vwind = x(14,:);wwind = x(15,:);
elseif num  == 9
    la_x = 0; la_y = 0; la_z = 0;
    la_p = 0; la_q = 0; la_r = 0;
elseif num  == 6
    la_x = 0; la_y = 0; la_z = 0;
    la_p = 0; la_q = 0; la_r = 0;
end
% input
Axm = uk(1,:); Aym = uk(2,:); Azm = uk(3,:);
pm = uk(4,:); qm = uk(5,:); rm = uk(6,:); 

% lameda_p=x(10);lameda_q=x(11);lameda_r=x(12);



F = ... 
[ 0, 0, 0,   (Aym - la_y)*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) + (Azm - la_z)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)), cos(phi)*cos(psi)*cos(theta)*(Azm - la_z) - cos(psi)*sin(theta)*(Axm - la_x) + cos(psi)*cos(theta)*sin(phi)*(Aym - la_y), (Azm - la_z)*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) - (Aym - la_y)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)) - cos(theta)*sin(psi)*(Axm - la_x);
  0, 0, 0, - (Aym - la_y)*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) - (Azm - la_z)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)), cos(phi)*cos(theta)*sin(psi)*(Azm - la_z) - sin(psi)*sin(theta)*(Axm - la_x) + cos(theta)*sin(phi)*sin(psi)*(Aym - la_y), (Azm - la_z)*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) - (Aym - la_y)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) + cos(psi)*cos(theta)*(Axm - la_x);
  0, 0, 0,                                                                 cos(phi)*cos(theta)*(Aym - la_y) - cos(theta)*sin(phi)*(Azm - la_z),                          - cos(theta)*(Axm - la_x) - cos(phi)*sin(theta)*(Azm - la_z) - sin(phi)*sin(theta)*(Aym - la_y),                                                                                                                                                                    0;
  0, 0, 0,                                                                   sin(phi)*tan(theta)*(la_r - rm) - cos(phi)*tan(theta)*(la_q - qm),                                      - cos(phi)*(la_r - rm)*(tan(theta)^2 + 1) - sin(phi)*(la_q - qm)*(tan(theta)^2 + 1),                                                                                                                                                                    0;
  0, 0, 0,                                                                                         cos(phi)*(la_r - rm) + sin(phi)*(la_q - qm),                                                                                                                        0,                                                                                                                                                                    0;
  0, 0, 0,                                                               (sin(phi)*(la_r - rm))/cos(theta) - (cos(phi)*(la_q - qm))/cos(theta),                        - (cos(phi)*sin(theta)*(la_r - rm))/cos(theta)^2 - (sin(phi)*sin(theta)*(la_q - qm))/cos(theta)^2,                                                                                                                                                                    0;];
 




