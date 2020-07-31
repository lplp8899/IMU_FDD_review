function F=F_lin(x,uk,gk)

% uGS_b=x(1);vGS_b=x(2);wGS_b=x(3);
phi = x(4); theta = x(5); psi = x(6);
% la_x=x(7);la_y=x(8);la_z=x(9);
[num,~]=size(x);
% if num==12
% %     la_x=x(7,:);la_y=x(8,:);la_z=x(9,:);
%     la_p=x(10,:);la_q=x(11,:);la_r=x(12,:);
% elseif num==9
%     la_V=x(7,:);la_a=x(8,:);la_b=x(9,:);
% %     la_x=0; la_y=0; la_z=0;
%     la_p=0; la_q=0; la_r=0;
% elseif num==6

    la_x=0; la_y=0; la_z=0;
    la_p=0; la_q=0; la_r=0;
% end
Axm = uk(1,:); Aym = uk(2,:); Azm = uk(3,:);
pm = uk(4,:); qm = uk(5,:); rm = uk(6,:); 
g = gk;

F = ... 
[ 0, 0, 0,   (Aym - la_y)*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) + (Azm - la_z)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)), cos(phi)*cos(psi)*cos(theta)*(Azm - la_z) - cos(psi)*sin(theta)*(Axm - la_x) + cos(psi)*cos(theta)*sin(phi)*(Aym - la_y), (Azm - la_z)*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) - (Aym - la_y)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)) - cos(theta)*sin(psi)*(Axm - la_x);
  0, 0, 0, - (Aym - la_y)*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) - (Azm - la_z)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)), cos(phi)*cos(theta)*sin(psi)*(Azm - la_z) - sin(psi)*sin(theta)*(Axm - la_x) + cos(theta)*sin(phi)*sin(psi)*(Aym - la_y), (Azm - la_z)*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) - (Aym - la_y)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) + cos(psi)*cos(theta)*(Axm - la_x);
  0, 0, 0,                                                                 cos(phi)*cos(theta)*(Aym - la_y) - cos(theta)*sin(phi)*(Azm - la_z),                          - cos(theta)*(Axm - la_x) - cos(phi)*sin(theta)*(Azm - la_z) - sin(phi)*sin(theta)*(Aym - la_y),                                                                                                                                                                    0;
  0, 0, 0,                                                                   sin(phi)*tan(theta)*(la_r - rm) - cos(phi)*tan(theta)*(la_q - qm),                                      - cos(phi)*(la_r - rm)*(tan(theta)^2 + 1) - sin(phi)*(la_q - qm)*(tan(theta)^2 + 1),                                                                                                                                                                    0;
  0, 0, 0,                                                                                         cos(phi)*(la_r - rm) + sin(phi)*(la_q - qm),                                                                                                                        0,                                                                                                                                                                    0;
  0, 0, 0,                                                               (sin(phi)*(la_r - rm))/cos(theta) - (cos(phi)*(la_q - qm))/cos(theta),                        - (cos(phi)*sin(theta)*(la_r - rm))/cos(theta)^2 - (sin(phi)*sin(theta)*(la_q - qm))/cos(theta)^2,                                                                                                                                                                    0;];


































