function E_n=E_lin(x)

%-------------------------  NOTE ----------------------------
% When adding to the model, it should be positive.
% If desiring the correct kinematic model, then it is negative.
% However, if use the positive, means the fault estimation is not used
% by the filter.
% use the SMO is different. It uses veq to compensate the faults.
%-------------------------  NOTE ----------------------------


% xGS=x(1);yGS=x(2);zGS=x(3);
% uGS_b=x(1);vGS_b=x(2);wGS_b=x(3);
phi = x(4); theta = x(5); psi = x(6);
% lameda_x=x(10);lameda_y=x(11);lameda_z=x(12);
% lameda_p=x(13);lameda_q=x(14);lameda_r=x(15);
% u_w=x(16);v_w=x(17);w_w=x(18);

[num,den]=size(x);
E_n = ...
[ -cos(psi)*cos(theta),   cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta), - sin(phi)*sin(psi) - cos(phi)*cos(psi)*sin(theta),  0,                    0,                    0;
  -cos(theta)*sin(psi), - cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta),   cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta),  0,                    0,                    0;
            sin(theta),                               -cos(theta)*sin(phi),                               -cos(phi)*cos(theta),  0,                    0,                    0;
                     0,                                                  0,                                                  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta);
                     0,                                                  0,                                                  0,  0,            -cos(phi),             sin(phi);
                     0,                                                  0,                                                  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta);];
                 
                 
% E_n = - E_n;

                 
    