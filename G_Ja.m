function G_n=G_Ja(x)

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


G_n = ...
[ -cos(psi)*cos(theta),   cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta), - sin(phi)*sin(psi) - cos(phi)*cos(psi)*sin(theta),  0,                    0,                    0;
  -cos(theta)*sin(psi), - cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta),   cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta),  0,                    0,                    0;
            sin(theta),                               -cos(theta)*sin(phi),                               -cos(phi)*cos(theta),  0,                    0,                    0;
                     0,                                                  0,                                                  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta);
                     0,                                                  0,                                                  0,  0,            -cos(phi),             sin(phi);
                     0,                                                  0,                                                  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta); ];
    